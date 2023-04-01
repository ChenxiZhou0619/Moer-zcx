#include "PathFilter.h"
#include <FunctionLayer/Material/Material.h>
#include <nanoflann/nanoflann.hpp>
#include <tbb/tbb.h>
void PathFilter::render(const Camera &camera, const Scene &scene,
                        std::shared_ptr<Sampler> sampler, int spp) const {
  int width = camera.film->size[0], height = camera.film->size[1], finished = 0;

  SVGroup group(width * height * spp);

  //* Sample the path first
  tbb::parallel_for(
      tbb::blocked_range2d<size_t>(0, width, 0, height),
      [&](const tbb::blocked_range2d<size_t> &r) {
        for (int row = r.rows().begin(); row != r.rows().end(); ++row) {
          for (int col = r.cols().begin(); col != r.cols().end(); ++col) {
            Vector2f NDC{(float)row / width, (float)col / height};
            for (int i = 0; i < spp; ++i) {
              ShadingVertex sv;
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
        }
      });
  printProgress(1);

  //* Filter the SVGroup

  using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
      nanoflann::L2_Simple_Adaptor<float, SVGroup>, SVGroup, 3 /* dim */
      >;

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

            std::vector<nanoflann::ResultItem<uint32_t, float>> ret_matches;

            size_t size =
                index.radiusSearch(&query_point[0], query_radius, ret_matches);

            //* Interpolate
            Spectrum Li_refine(.0f);
            float Li_weight = .0f;
            for (int i = 0; i < size; ++i) {
              size_t idx = ret_matches[i].first;
              const auto &xj = group.svs[idx];
              if (float weight = dot(sv.normal, xj.normal); weight > .8f) {
                Li_refine += xj.Li * weight;
                Li_weight += weight;
              }
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

void PathFilter::samplePath(Ray ray, const Scene &scene,
                            std::shared_ptr<Sampler> sampler,
                            ShadingVertex *sv) const {
  Spectrum Li(.0f), Le(.0f), beta(1.f);
  int depth = 0;

  while (true) {
    auto si = scene.rayIntersect(ray);
    Vector3f wo = -ray.direction;

    if (depth == 0 /* si respond to a shading vertex */) {
      if (!si /* This shading point cann't be filter */) {
        sv->valid = false;
        sv->weight = Spectrum(1.f);
        for (auto light : scene.infiniteLights) {
          Le += light->evaluateEmission(ray);
        }
        break;
      } else if (auto light = si->shape->light; light) {
        Le += light->evaluateEmission(*si, wo);
      }
      sv->valid = true;
      sv->position = si->position;
      sv->normal = si->normal;
      sv->weight = Spectrum(1.f);
      sv->distSqr = si->distance * si->distance;
    }

    if (++depth > maxDepth)
      break;

    if (!si)
      break;

    //* Compute incident radiance for shading vertex vs
    auto material = si->shape->material;
    auto bsdf = material->computeBSDF(*si);

    //* Next Event Estimation

    //* Sample infinite lights
    for (auto light : scene.infiniteLights) {
      auto result = light->sample(*si, sampler->next2D());
      Ray shadowRay{si->position, result.direction, 1e-4f, result.distance};
      if (auto occlude = scene.rayIntersect(shadowRay); !occlude) {
        Vector3f wi = shadowRay.direction;
        Spectrum f = bsdf->f(wo, wi);
        float pdf = convertPDF(result, *si),
              misw =
                  result.isDelta ? 1.f : powerHeuristic(pdf, bsdf->pdf(wo, wi));
        if (!f.isZero() && pdf != .0f) {
          Li += beta * f * result.energy * misw / pdf;
        }
      }
    }

    //* Sample light in scene
    float pdfLight = .0f;
    if (auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
        light && pdfLight != .0f) {
      auto result = light->sample(*si, sampler->next2D());
      Ray shadowRay{si->position, result.direction, 1e-4f, result.distance};
      if (auto occlude = scene.rayIntersect(shadowRay); !occlude) {
        Vector3f wi = shadowRay.direction;
        Spectrum f = bsdf->f(wo, wi);
        float pdf = convertPDF(result, *si) * pdfLight,
              misw =
                  result.isDelta ? 1.f : powerHeuristic(pdf, bsdf->pdf(wo, wi));
        if (!f.isZero() && pdf != .0f) {
          Li += beta * f * result.energy * misw / pdf;
        }
      }
    }

    //* Sample BSDF
    {
      auto result = bsdf->sample(wo, sampler->next2D());
      Ray shadowRay{si->position, result.wi, 1e-4f, FLT_MAX};
      auto Le_si = scene.rayIntersect(shadowRay);
      if (!Le_si) {
        for (auto light : scene.infiniteLights) {
          float pdf = result.pdf,
                misw = powerHeuristic(pdf, light->pdf(shadowRay));
          Li +=
              beta * result.weight * light->evaluateEmission(shadowRay) * misw;
        }
      } else if (auto light = Le_si->shape->light; light) {
        float pdf = result.pdf,
              pdfLe = scene.pdf(light) * light->pdf(shadowRay, *Le_si),
              misw = powerHeuristic(pdf, pdfLe);
        Li += beta * result.weight *
              light->evaluateEmission(*Le_si, -shadowRay.direction) * misw;
      }
    }

    //* Spwan the ray
    auto result = bsdf->sample(wo, sampler->next2D());
    ray = Ray{si->position, result.wi, 1e-4f, FLT_MAX};
    beta *= result.weight;
  }

  sv->Li = Li;
  sv->Le = Le;
}

REGISTER_CLASS(PathFilter, "pathFilter")