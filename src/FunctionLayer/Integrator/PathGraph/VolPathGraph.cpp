#include "VolPathGraph.h"
#include <FunctionLayer/Material/Material.h>
#include <array>
#include <nanoflann/nanoflann.hpp>
#include <stack>
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

  auto propagation = [&]() {
    //* Propagation once
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, width * height),
        [&](const tbb::blocked_range<size_t> &r) {
          for (int i = r.begin(); i != r.end(); ++i) {
            int index = i;
            //* index = x + y * width
            int y = i / width, x = i - y * width;
            for (int j = 0; j < spp; ++j) {
              if (int idx = fspGroup[i * spp + j]; idx != -1) {
                VolumeShadingPoint *sp = &group.points[idx];

                std::stack<VolumeShadingPoint *> path;

                while (sp != nullptr) {
                  path.push(sp);
                  sp = sp->indirect.xj;
                }

                while (!path.empty()) {
                  VolumeShadingPoint *vs = path.top();
                  path.pop();

                  if (vs->indirect.xj)
                    vs->indirect.li = vs->indirect.xj->lo;

                  //* Compute lo of vs
                  Spectrum lo{.0f};

                  // nee
                  lo += vs->nee.weight * vs->nee.fp * vs->nee.li / vs->nee.pdf;
                  // phs
                  lo += vs->phs.weight * vs->phs.fp * vs->phs.li / vs->phs.pdf;
                  // indirect
                  if (vs->indirect.xj)
                    lo += vs->indirect.fp * vs->indirect.li / vs->indirect.pdf;

                  vs->lo = lo;
                }
              }
            }
          }
        });
  };

  using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<float, VSPGroup>, VSPGroup, 3 /* dim */
      >;
  my_kd_tree_t index(3, group, 32);

  finished = 0;
  std::cout << "\n";
  // Compute neighborhoods
  tbb::parallel_for(
      tbb::blocked_range<size_t>(0, group.points.size()),
      [&](const tbb::blocked_range<size_t> &r) {
        for (int i = r.begin(); i != r.end(); ++i) {
          auto &p = group.points[i];

          float query_point[3] = {p.position[0], p.position[1], p.position[2]};
          std::vector<uint32_t> ret_index(knn);
          std::vector<float> out_dist_sqr(knn);

          size_t num = index.knnSearch(query_point, knn, &ret_index[0],
                                       &out_dist_sqr[0]);
          for (int j = 0; j < num; ++j) {
            p.neighbors.emplace_back(&group.points[ret_index[j]]);
          }

          //* Compute mis weight for each neighbor
          for (int j = 0; j < num; ++j) {
            float nee_rho = .0f;
            float phs_rho = .0f;
            float ind_rho = .0f;

            Vector3f nee_wi = p.neighbors[j]->nee.wi,
                     phs_wi = p.neighbors[j]->phs.wi,
                     ind_wi = p.neighbors[j]->indirect.wi;

            for (int l = 0; l < num; ++l) {
              VolumeShadingPoint *neighbor = p.neighbors[l];
              nee_rho += neighbor->phase->pdf(neighbor->wo, nee_wi);
              nee_rho += scene.infiniteLights[0]->pdf(neighbor->shadowRay_nee);

              phs_rho += neighbor->phase->pdf(neighbor->wo, phs_wi);
              phs_rho += scene.infiniteLights[0]->pdf(neighbor->shadowRay_phs);

              if (neighbor->indirect.xj) {
                ind_rho +=
                    neighbor->phase->pdf(neighbor->wo, neighbor->indirect.wi);
              }
            }

            p.nee_rho.emplace_back(1.f / nee_rho);
            p.phs_rho.emplace_back(1.f / phs_rho);

            if (ind_rho > .0f)
              p.ind_rho.emplace_back(1.f / ind_rho);
            else {
              p.ind_rho.emplace_back(.0f);
            }
          }

          ++finished;
          if (finished % 20 == 0) {
            printProgress((float)finished / group.points.size());
          }
        }
      });

  for (int i = 0; i < iterations; ++i) {

    propagation();

    tbb::parallel_for(tbb::blocked_range<size_t>(0, group.points.size()),
                      [&](const tbb::blocked_range<size_t> &r) {
                        for (int i = r.begin(); i != r.end(); ++i) {
                          Spectrum Lo{.0f};
                          VolumeShadingPoint *p = &group.points[i];
                          Vector3f wo = p->wo;
                          auto phase = p->phase;
                          //* Refine p->lo
                          for (int j = 0; j < p->neighbors.size(); ++j) {
                            VolumeShadingPoint *neighbor = p->neighbors[j];
                            Lo += phase->f(wo, neighbor->nee.wi) *
                                  p->nee_rho[j] * neighbor->nee.li;
                            Lo += phase->f(wo, neighbor->phs.wi) *
                                  p->phs_rho[j] * neighbor->phs.li;
                          }
                          //* Refine lo through indirect
                          if (p->indirect.xj) {
                            for (int j = 0; j < p->neighbors.size(); ++j) {
                              VolumeShadingPoint *neighbor = p->neighbors[j];
                              if (p->ind_rho[j] > .0f) {
                                Lo += phase->f(wo, neighbor->indirect.wi) *
                                      neighbor->indirect.li * p->ind_rho[j];
                              }
                            }
                          }
                          p->lo = Lo;
                        }
                      });
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
                            Spectrum L = sp.lo;
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
          vsp.nee.fp = p;
          vsp.nee.pdf = pdf;
          vsp.nee.li = Tr * result.energy;
          vsp.shadowRay_nee = shadwoRay;
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
              vsp.phs.pdf = result.pdf;
              vsp.phs.fp = phase->f(wo, wi);
              vsp.phs.li = Tr * energy;

              vsp.shadowRay_phs = shadowRay;
            }
          }
        }

        PhaseSampleResult result = phase->sample(wo, sampler->next2D());
        ray = Ray{mi.position, result.wi, 1e-4f, FLT_MAX};
        ray.medium = medium;

        vsp.indirect.fp = phase->f(wo, result.wi);
        vsp.indirect.pdf = phase->pdf(wo, result.wi);
        vsp.indirect.wi = result.wi;

        int idx = group->emplace_back(vsp);
        if (depth == 1 /** First shading point */) {
          *firstShadingPointIdx = idx;
        }
        if (prevShadingPointIdx != -1) {
          //* Set the prev point to current point if prev exists
          group->points[prevShadingPointIdx].indirect.xj = &group->points[idx];
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