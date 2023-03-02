// PathIntegrator with mis and volume rendering

#include "PathIntegrator.h"
#include <FunctionLayer/Material/Material.h>

float powerHeuristic(float pdfA, float pdfB) {
  pdfA *= pdfA;
  pdfB *= pdfB;
  return pdfA / (pdfA + pdfB);
}

Spectrum PathIntegrator::li(const Ray &_ray, const Scene &scene,
                            std::shared_ptr<Sampler> sampler) const {
  Spectrum spectrum(.0f), beta(1.f);
  Ray ray(_ray);
  int step = 0;

  float pdfPrev;
  bool isDeltaPrev = true;

  do {
    auto itsOpt = scene.rayIntersect(ray);

    if (!itsOpt.has_value()) {
      for (auto light : scene.infiniteLights) {
        float pdf = light->pdf(ray);
        float misw = isDeltaPrev ? 1.f : powerHeuristic(pdfPrev, pdf);
        spectrum += beta * misw * light->evaluateEmission(ray);
      }
      break; // Break due to escaping the scene
    }

    Intersection its = itsOpt.value();
    if (auto light = its.shape->light; light) {
      float pdf = light->pdf(ray, its);
      pdf *= scene.pdf(light);
      float misw = isDeltaPrev ? 1.f : powerHeuristic(pdfPrev, pdf);
      spectrum += beta * misw * light->evaluateEmission(its, -ray.direction);
    }

    // TODO RR to determine whether to stop the ray
    if (step > 3)
      break;

    computeRayDifferentials(&its, ray);
    auto bsdf = its.shape->material->computeBSDF(its);

    // If hit the empty surface
    if (!bsdf) {
      ray = Ray(its.position, ray.direction);
      continue;
    }

    // Sample the InfiniteLights
    for (auto light : scene.infiniteLights) {
      auto lightSampleResult = light->sample(its, sampler->next2D());
      Ray shadowRay(its.position, lightSampleResult.direction, 1e-4f,
                    lightSampleResult.distance);
      if (auto occlude = scene.rayIntersect(shadowRay); !occlude.has_value()) {
        Spectrum f = bsdf->f(-ray.direction, shadowRay.direction);
        float pdf = convertPDF(lightSampleResult, its);
        float misw = lightSampleResult.isDelta
                         ? 1.f
                         : powerHeuristic(pdf, bsdf->pdf(-ray.direction,
                                                         shadowRay.direction));
        if (!f.isZero() && pdf != .0f)
          spectrum += beta * misw * lightSampleResult.energy * f / pdf;
      }
    }

    // Sample the lights in scene
    float pdfLight = .0f;
    auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
    if (light && pdfLight != .0f) {
      auto lightSampleResult = light->sample(its, sampler->next2D());
      Ray shadowRay(its.position, lightSampleResult.direction, 1e-4f,
                    lightSampleResult.distance);
      if (auto occlude = scene.rayIntersect(shadowRay); !occlude.has_value()) {
        Spectrum f = bsdf->f(-ray.direction, shadowRay.direction);
        lightSampleResult.pdf *= pdfLight;
        float pdf = convertPDF(lightSampleResult, its);
        float misw = lightSampleResult.isDelta
                         ? 1.f
                         : powerHeuristic(pdf, bsdf->pdf(-ray.direction,
                                                         shadowRay.direction));
        if (!f.isZero() && pdf != .0f)
          spectrum += beta * misw * lightSampleResult.energy * f / pdf;
      }
    }

    // Sample the bsdf to spwan the ray
    auto bsdfSampleResult = bsdf->sample(-ray.direction, sampler->next2D());
    beta *= bsdfSampleResult.weight;
    if (beta.isZero())
      break;
    ray = Ray(its.position, bsdfSampleResult.wi, 1e-4f, FLT_MAX);

    isDeltaPrev = bsdfSampleResult.type == BSDFType::Specular;
    pdfPrev = bsdfSampleResult.pdf;
    ++step;

  } while (1);

  return spectrum;
}

REGISTER_CLASS(PathIntegrator, "path")