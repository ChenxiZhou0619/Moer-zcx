#include "VolPathFilter.h"
#include <FunctionLayer/Material/Material.h>
#include <nanoflann/nanoflann.hpp>
#include <tbb/tbb.h>

int SampleMode1(float sigma_a, float sigma_s, float sigma_n, float uMode) {
  uMode *= (sigma_a + sigma_s + sigma_n);
  if (uMode < sigma_a)
    return 0; // Absorb
  else if (uMode < sigma_a + sigma_s)
    return 1; // Scatter
  else
    return 2; // Null-collision
}

Spectrum VolPathFilter::Transmittance(
    const Scene &scene, Ray ray,
    std::function<Spectrum(const Medium *, Ray, float)> tr_tracker) const {
  Spectrum Tr(1.f);

  while (true) {
    const Medium *medium = ray.medium;

    auto si = scene.rayIntersect(ray);
    float dt = si ? si->distance : ray.tFar;

    if (medium)
      Tr *= tr_tracker(medium, ray, dt);

    if (si) {
      auto bsdf = si->shape->material->computeBSDF(*si);
      if (bsdf) {
        //* Occlude
        Tr = Spectrum(.0f);
        break;
      } else {
        //* Skip this surface
        ray = Ray{si->position, ray.direction, 1e-4f, ray.tFar - dt};
        setRayMedium(ray.direction, si->normal, si->shape->material, &ray);
        continue;
      }
    }
    break;
  }

  return Tr;
}

std::optional<Intersection> VolPathFilter::TransmittanceRayIntersect(
    const Scene &scene, Ray ray, Spectrum *Tr,
    std::function<Spectrum(const Medium *, Ray, float)> tr_tracker) const {
  float distance = .0f;

  while (true) {
    const Medium *medium = ray.medium;
    auto si = scene.rayIntersect(ray);

    float dt = si ? si->distance : ray.tFar;

    if (medium)
      *Tr *= tr_tracker(medium, ray, dt);

    if (!si)
      return std::nullopt;

    auto material = si->shape->material;
    auto bsdf = material->computeBSDF(*si);
    distance += dt;

    if (bsdf) {
      //* Return the surface intersection
      si->distance = distance;
      return si;
    } else {
      //* Skip this surface
      ray = Ray{si->position, ray.direction, 1e-4f, ray.tFar - dt};
      setRayMedium(ray.direction, si->normal, material, &ray);
      continue;
    }
  }
}

void VolPathFilter::render(const Camera &camera, const Scene &scene,
                           std::shared_ptr<Sampler> sampler, int spp) const {
  int width = camera.film->size[0], height = camera.film->size[1];
  int finished = 0;

  VolSVGroup group(width * height * spp);

  std::cout << "VolPathFilter render!\n";

  //* Sample path first
  tbb::parallel_for(
      tbb::blocked_range2d<size_t>(0, width, 0, height),
      [&](const tbb::blocked_range2d<size_t> &r) {
        for (int row = r.rows().begin(); row != r.rows().end(); ++row)
          for (int col = r.cols().begin(); col != r.cols().end(); ++col) {
            Vector2f NDC{(float)row / width, (float)col / height};
            for (int i = 0; i < spp; ++i) {
              VolShadingVertex sv;

              sv.pixelLoc = Vector2i{row, col};
              sv.pixelSA = camera.pixelSolidAngle(NDC);

              Ray ray = camera.sampleRayDifferentials(
                  CameraSample{sampler->next2D()}, NDC);

              samplePath(ray, scene, sampler, &sv);

              if (sv.valid)
                group.emplace_back(sv);
              else {
                Spectrum L = sv.weight * (sv.Le + sv.Li);
                camera.film->deposit({row, col}, L, 1.f);
              }
            }
            ++finished;

            if (finished % 20 == 0) {
              printProgress((float)finished / (height * width));
            }
          }
      });
  printProgress(1.f);

  using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<float, VolSVGroup>, VolSVGroup, 3>;

  my_kd_tree_t index(3, group, 32);

  const float n = fm::pow<float>(group.svs.size(), .25f);

  for (int itr = 0; itr < iteration; ++itr) {
    //* Compute refine
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, group.svs.size()),
        [&](const tbb::blocked_range<size_t> &r) {
          for (int i = r.begin(); i != r.end(); ++i) {
            auto &sv = group.svs[i];
            float query_point[3];

            query_point[0] = sv.position[0];
            query_point[1] = sv.position[1];
            query_point[2] = sv.position[2];

            float query_radius = fm::sqrt(sv.distSqr * sv.pixelSA * INV_PI) / n;

            query_radius *= k;

            std::vector<nanoflann::ResultItem<uint32_t, float>> ret_matches;
            size_t size =
                index.radiusSearch(&query_point[0], query_radius, ret_matches);

            //* Interpolate
            Spectrum Li_refine(.0f);
            float Li_weight = .0f;
            for (int i = 0; i < size; ++i) {
              size_t idx = ret_matches[i].first;
              const auto &xj = group.svs[idx];

              //              float weight = 1.f; // TODO set the weight here

              float weight = sv.sigma_t / xj.sigma_t;

              if (weight > 1.f)
                weight = 1.f / weight;

              float dist = (xj.position - sv.position).length();
              weight *= std::exp(-dist / query_radius);

              Li_refine += xj.Li * weight;
              Li_weight += weight;
            }

            sv.Li_refine = Li_refine / Li_weight;
          }
        });

    //* Update refine
    tbb::parallel_for(tbb::blocked_range<size_t>(0, group.svs.size()),
                      [&](const tbb::blocked_range<size_t> &r) {
                        for (int i = r.begin(); i != r.end(); ++i) {
                          auto &sv = group.svs[i];
                          sv.Li = sv.Li_refine;
                        }
                      });
  }

  //* Reconstruct the image
  for (const auto &sv : group.svs) {
    Spectrum li = sv.weight * (sv.Le + sv.Li);
    camera.film->deposit(sv.pixelLoc, li, 1.f);
  }
}

void VolPathFilter::samplePath(Ray ray, const Scene &scene,
                               std::shared_ptr<Sampler> sampler,
                               VolShadingVertex *sv) const {
  Spectrum Li(.0f), Le(.0f), beta(1.f);
  float firstLen = .0f;
  int depth = 0;

  auto tr_tracker_regular = [](const Medium *medium, Ray ray, float tmax) {
    return medium->Transmittance_RatioTracking(ray, tmax);
  };

  while (true) {
    auto si = scene.rayIntersect(ray);
    Vector3f wo = -ray.direction;

    const Medium *medium = ray.medium;
    MediumIntersection mi;
    bool mediumInteraction = false;

    if (medium) {
      Spectrum Tr;
      float pdf, tmax = si ? si->distance : FLT_MAX;
      mediumInteraction = medium->Sample_MajorantTracking(
          ray, tmax, sampler->next2D(), &mi, &Tr, &pdf);
    }

    if (mediumInteraction) {

      //* Handle medium interaction
      if (++depth > maxDepth)
        break;

      firstLen += mi.distance;

      auto phase = mi.mp.phase;

      float sigma_a = mi.mp.sigma_a[0], sigma_s = mi.mp.sigma_s[0],
            sigma_n = std::max(.0f, mi.mp.sigma_maj[0] - sigma_a - sigma_s);

      int mode = SampleMode1(sigma_a, sigma_s, sigma_n, sampler->next1D());

      if (mode == 0 /** Absorb */) {
        // Just terminate
        // TODO
        break;
      } else if (mode == 1 /** Scatter */) {
        // Record the first diffuse scatter point

        if (depth == 1 /** mi respond to a shading vertex */) {
          sv->valid = true;
          sv->weight = Spectrum(1.f);

          sv->position = mi.position;
          sv->distSqr = firstLen * firstLen;

          sv->sigma_t = mi.mp.sigma_a[0] + mi.mp.sigma_s[0];
        }

        //* 1. Sample Le
        for (auto light : scene.infiniteLights) {
          LightSampleResult result = light->sample(mi, sampler->next2D());

          Ray shadowRay(mi.position, result.direction, 1e-4f, result.distance);
          shadowRay.medium = medium;

          Vector3f wi = shadowRay.direction;
          Spectrum Tr = Transmittance(scene, shadowRay, tr_tracker_regular),
                   p = phase->f(wo, wi);
          float pdf = convertPDF(result, mi),
                misw = result.isDelta ? 1.f
                                      : powerHeuristic(pdf, phase->pdf(wo, wi));

          if (!Tr.isZero() && !p.isZero()) {
            Li += beta * p * Tr * result.energy * misw / pdf;
          }
        }

        //* 2. Sample fp
        {
          PhaseSampleResult result = phase->sample(wo, sampler->next2D());

          Ray shadowRay{mi.position, result.wi, 1e-4f, FLT_MAX};
          shadowRay.medium = medium;
          Vector3f wi = shadowRay.direction;

          Spectrum Tr(1.f), p = phase->f(wo, wi);
          auto Le_si = TransmittanceRayIntersect(scene, shadowRay, &Tr,
                                                 tr_tracker_regular);
          if (!Le_si) {
            //* Hit environment
            for (auto light : scene.infiniteLights) {
              float pdf = result.pdf,
                    misw = result.isDelta
                               ? 1.f
                               : powerHeuristic(pdf, light->pdf(shadowRay));
              Spectrum energy = light->evaluateEmission(shadowRay);
              Li += beta * result.weight * Tr * energy * misw;
            }
          }
        }

        PhaseSampleResult result = phase->sample(wo, sampler->next2D());

        ray = Ray{mi.position, result.wi, 1e-4f, FLT_MAX};
        ray.medium = medium;
        beta *= result.weight;

      } else /** Null-collision */ {
        --depth;

        ray = Ray{mi.position, ray.direction, 1e-4f, FLT_MAX};
        ray.medium = medium;
      }
    } else {
      //* Handle surface interaction
      //* Cuz the test scene not contains surface (except empty)
      //* So, just skip the empty and terminate(error) else

      if (depth == 0) {
        if (!si) {
          for (auto light : scene.infiniteLights)
            Le += beta * light->evaluateEmission(ray);
        } else if (auto light = si->shape->light; light) {
          Le += beta * light->evaluateEmission(*si, wo);
        }
      }
      if (!si) {
        sv->valid = false;
        break;
      }

      auto material = si->shape->material;
      auto bsdf = material->computeBSDF(*si);

      firstLen += si->distance;

      if (!bsdf) {
        // Just skip the surface
        ray = Ray(si->position, ray.direction, 1e-4f, FLT_MAX);
        setRayMedium(ray.direction, si->normal, material, &ray);
        continue;
      } else {
        // !Error, there shouldn't be any not empty surface
        std::cout << "Error!\n";
        std::exit(1);
      }
    }
  }
  sv->weight = Spectrum(1.f);
  sv->Li = Li;
  sv->Le = Le;
}

void VolPathFilter::setRayMedium(Vector3f direction, Vector3f normal,
                                 std::shared_ptr<Material> material,
                                 Ray *ray) const {
  bool towardsInner = dot(direction, normal) < .0f;
  ray->medium = towardsInner ? material->getMedium() : nullptr;
}

REGISTER_CLASS(VolPathFilter, "volPathFilter")