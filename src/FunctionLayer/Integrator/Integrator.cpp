#include "Integrator.h"
#include <tbb/tbb.h>

#define MULTITHREAD

Spectrum Integrator::sampleDirect(const Scene &scene, const Intersection &its,
                                  Vector3f wo, std::shared_ptr<BSDF> bsdf,
                                  std::shared_ptr<Sampler> sampler) const {
  Spectrum Ld{.0f};

  // Sample Le
  auto sampleLd = [&](std::shared_ptr<Light> light, Vector2f u2,
                      float pdfLight = 1.f) {
    auto result = light->sample(its, u2);
    Vector3f wi = result.direction;
    Ray shadowRay{its.position, wi, 1e-4f, result.distance};

    if (auto occlusion = scene.rayIntersect(shadowRay); !occlusion) {
      Spectrum f = bsdf->f(wo, wi);
      float pdf = convertPDF(result, its) * pdfLight;
      float misw =
          result.isDelta ? 1.f : powerHeuristic(pdf, bsdf->pdf(wo, wi));
      if (!f.isZero() && pdf != .0f)
        return misw * result.energy * f / pdf;
    }
    return Spectrum{.0f};
  };
  /** ----- Infinite Lights ----- */
  for (auto light : scene.infiniteLights) {
    Ld += sampleLd(light, sampler->next2D(), 1.f);
  }
  /** ----- Surface Lights ----- */
  float pdfLight = .0f;
  auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
  if (light && pdfLight != .0f)
    sampleLd(light, sampler->next2D(), pdfLight);

  // Sample bsdf
  auto result = bsdf->sample(wo, sampler->next2D());
  Ray shadowRay{its.position, result.wi, 1e-4f, FLT_MAX};
  auto occlusion = scene.rayIntersect(shadowRay);

  if (!occlusion) {
    for (auto light : scene.infiniteLights) {
      float pdfLe = light->pdf(shadowRay);
      float misw = (result.type == BSDFType::Specular)
                       ? 1.f
                       : powerHeuristic(result.pdf, pdfLe);
      Ld += result.weight * misw * light->evaluateEmission(shadowRay);
    }
  } else if (auto light = occlusion->shape->light; light) {
    float pdfLight = scene.pdf(light);
    float pdfLe = pdfLight * light->pdf(shadowRay, *occlusion);
    float misw = (result.type == BSDFType::Specular)
                     ? 1.f
                     : powerHeuristic(result.pdf, pdfLe);
    Ld +=
        result.weight * misw * light->evaluateEmission(*occlusion, -result.wi);
  }

  return Ld;
}

void PixelIntegrator::render(const Camera &camera, const Scene &scene,
                             std::shared_ptr<Sampler> sampler, int spp) const {
  int width = camera.film->size[0], height = camera.film->size[1], finished = 0;
#ifdef MULTITHREAD
  tbb::parallel_for(
      tbb::blocked_range2d<size_t>(0, width, 0, height),
      [&](const tbb::blocked_range2d<size_t> &r) {
        for (int row = r.rows().begin(); row != r.rows().end(); ++row)
          for (int col = r.cols().begin(); col != r.cols().end(); ++col) {
            Vector2f NDC{(float)row / width, (float)col / height};
            for (int i = 0; i < spp; ++i) {
              Ray ray = camera.sampleRayDifferentials(
                  CameraSample{sampler->next2D()}, NDC);
              Spectrum res = li(ray, scene, sampler);
              if (!res.hasNaN() && !res.hasInf())
                camera.film->deposit({row, col}, res, 1.f);
            }

            ++finished;
            if (finished % 20 == 0) {
              printProgress((float)finished / (height * width));
            }
          }
      });
#else
  for (int row = 0; row != width; ++row)
    for (int col = 0; col != height; ++col) {
      Vector2f NDC{(float)row / width, (float)col / height};
      for (int i = 0; i < spp; ++i) {
        Ray ray =
            camera.sampleRayDifferentials(CameraSample{sampler->next2D()}, NDC);
        camera.film->deposit({row, col}, li(ray, scene, sampler), 1.f);
      }

      ++finished;
      if (finished % 20 == 0) {
        printProgress((float)finished / (height * width));
      }
    }
#endif
  printProgress(1);
}