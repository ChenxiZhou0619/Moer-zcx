#include "VolPathGraph.h"
#include <FunctionLayer/Material/Material.h>
#include <array>
#include <tbb/tbb.h>

Spectrum VolPathGraph::Transmittance(
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

std::optional<Intersection> VolPathGraph::TransmittanceRayIntersect(
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

void VolPathGraph::setRayMedium(Vector3f direction, Vector3f normal,
                                std::shared_ptr<Material> material,
                                Ray *ray) const {
  bool towardsInner = dot(direction, normal) < .0f;
  ray->medium = towardsInner ? material->getMedium() : nullptr;
}

void VolPathGraph::render(const Camera &camera, const Scene &scene,
                          std::shared_ptr<Sampler> sampler, int spp) const {
  int width = camera.film->size[0], height = camera.film->size[1];

  VSPGroup group(width * height * spp * maxPathLength);

  std::vector<int> fspGroup(width * height * spp);
  std::fill(fspGroup.begin(), fspGroup.end(), -1);

  int finished = 0;

  //* Generate path graph
  tbb::parallel_for(
      tbb::blocked_range2d<size_t>(0, width, 0, height),
      [&](const tbb::blocked_range2d<size_t> &r) {
        for (int row = r.rows().begin(); row != r.rows().end(); ++row)
          for (int col = r.cols().begin(); col != r.cols().end(); ++col) {
            Vector2f NDC{(float)row / width, (float)col / height};
            for (int i = 0; i < spp; ++i) {
              Ray ray = camera.sampleRayDifferentials(
                  CameraSample{sampler->next2D()}, NDC);

              int firstShadingPointIdx =
                  -1; // -1 represents current path with no fsp
              Spectrum li = constructPath(ray, scene, sampler, &group,
                                          &firstShadingPointIdx);

              if (firstShadingPointIdx == -1) {
                camera.film->deposit({row, col}, li, 1);
              }

              fspGroup[(col * width + row) * spp + i] = firstShadingPointIdx;
            }
            ++finished;
            if (finished % 20 == 0) {
              printProgress((float)finished / (height * width));
            }
          }
      });
  printProgress(1);

  for (int i = 0; i < iterations; ++i) {
    //* Refine the direct incident estimation by filtering

    //* Refine the indirect incident estimation by propagation
  }

  //* Reconstruction the image
  tbb::parallel_for(tbb::blocked_range<size_t>(0, width * height),
                    [&](const tbb::blocked_range<size_t> &r) {
                      for (int i = r.begin(); i != r.end(); ++i) {
                        int index = i;
                        //* index = x + y * width
                        int y = i / width, x = i - y * width;
                        for (int j = 0; j < spp; ++j) {
                          if (int idx = fspGroup[i * spp + j]; idx != -1) {
                            auto sp = group.points[idx];
                            Spectrum L = sp.nee.li * sp.nee.weight +
                                         sp.phs.li * sp.phs.weight + sp.lj;
                            camera.film->deposit({x, y}, L, 1);
                          }
                        }
                      }
                    });
}

Spectrum VolPathGraph::constructPath(Ray ray, const Scene &scene,
                                     std::shared_ptr<Sampler> sampler,
                                     VSPGroup *group,
                                     int *firstShadingPointIdx) const {
  // When there is no shading point for current pixelLoc
  // Just return the li of the current ray

  Spectrum Li(.0f);
  int depth = 0;
  int prevShadingPointIdx = -1;

  auto tr_tracker = [](const Medium *medium, Ray ray, float tmax) {
    return medium->Transmittance_RatioTracking(ray, tmax);
  };

  auto sample_mode = [](float sigma_a, float sigma_s, float sigma_n,
                        float uMode) {
    uMode *= (sigma_a + sigma_s + sigma_n);
    if (uMode < sigma_a)
      return 0; // Absorb
    else if (uMode < sigma_a + sigma_s)
      return 1; // Scatter
    else
      return 2; // Null-collision
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
      if (++depth > maxPathLength)
        break;

      auto phase = mi.mp.phase;
      float sigma_a = mi.mp.sigma_a[0], sigma_s = mi.mp.sigma_s[0],
            sigma_n = std::max(.0f, mi.mp.sigma_maj[0] - sigma_a - sigma_s);
      int mode = sample_mode(sigma_a, sigma_s, sigma_n, sampler->next1D());

      VolumeShadingPoint vsp;

      if (mode == 0 /**Absorb  */) {
        // Just terminate
        break;
      } else if (mode == 1 /** Scatter */) {

        // Fill vsp
        vsp.wo = wo;
        vsp.phase = phase;
        vsp.position = mi.position;

        //* Sample Le direct
        for (auto light : scene.infiniteLights) {
          LightSampleResult result = light->sample(mi, sampler->next2D());

          Ray shadwoRay(mi.position, result.direction, 1e-4f, result.distance);
          shadwoRay.medium = medium;

          Vector3f wi = shadwoRay.direction;
          Spectrum Tr = Transmittance(scene, shadwoRay, tr_tracker),
                   p = phase->f(wo, wi);
          float pdf = convertPDF(result, mi),
                weight = result.isDelta
                             ? 1.f
                             : powerHeuristic(pdf, phase->pdf(wo, wi));

          vsp.nee.wi = wi;
          vsp.nee.weight = weight;
          vsp.nee.li = p * Tr * result.energy / pdf;
        }
        //* Sample phase direct
        {
          PhaseSampleResult result = phase->sample(wo, sampler->next2D());
          Ray shadowRay{mi.position, result.wi, 1e-4f, FLT_MAX};
          shadowRay.medium = medium;

          Vector3f wi = shadowRay.direction;
          Spectrum Tr(1.f), p = result.weight;
          auto Le_si =
              TransmittanceRayIntersect(scene, shadowRay, &Tr, tr_tracker);
          if (!Le_si) {
            //* Hit environment
            for (auto light : scene.infiniteLights) {
              float pdf = result.pdf,
                    weight = result.isDelta
                                 ? 1.f
                                 : powerHeuristic(pdf, light->pdf(shadowRay));
              Spectrum energy = light->evaluateEmission(shadowRay);
              vsp.phs.wi = wi;
              vsp.phs.weight = weight;
              vsp.phs.li = p * Tr * energy;
            }
          }
        }

        PhaseSampleResult result = phase->sample(wo, sampler->next2D());
        ray = Ray{mi.position, result.wi, 1e-4f, FLT_MAX};
        ray.medium = medium;

        vsp.wj = result.wi;
        vsp.weight_j = result.weight;

        int idx = group->emplace_back(vsp);
        if (depth == 1 /** First shading point */) {
          *firstShadingPointIdx = idx;
        }
        if (prevShadingPointIdx != -1) {
          //* Set the prev point to current point if prev exists
          group->points[prevShadingPointIdx].xj = &group->points[idx];
        }
        prevShadingPointIdx = idx;
      } else if (mode == 2 /** Null collision */) {
        --depth;

        ray = Ray{mi.position, ray.direction, 1e-4f, FLT_MAX};
        ray.medium = medium;
        continue;
      } else {
        std::cerr << "Error, shouldn't arrive here!\n";
      }
    } else {
      //* Surface interaction
      //* Only handle empty surface

      if (depth == 0) {
        if (!si) {
          for (auto light : scene.infiniteLights) {
            Li += light->evaluateEmission(ray);
          }
        }
      }

      if (!si) {
        break;
      }

      auto material = si->shape->material;
      auto bsdf = material->computeBSDF(*si);

      if (!bsdf) {
        ray = Ray(si->position, ray.direction, 1e-4f, FLT_MAX);
        setRayMedium(ray.direction, si->normal, material, &ray);
        continue;
      } else {
        // !Error, there shouldn't be any not empty surface
        std::cout << "Error! Hit not empty surface\n";
        std::exit(1);
      }
    }
  }

  return Li;
}

REGISTER_CLASS(VolPathGraph, "volPathGraph")