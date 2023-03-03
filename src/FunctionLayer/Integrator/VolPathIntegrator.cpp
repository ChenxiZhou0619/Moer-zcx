#include "VolPathIntegrator.h"
#include <FunctionLayer/Material/Material.h>
#include <FunctionLayer/Medium/Medium.h>

Spectrum VolPathIntegrator::Tr(const Scene &scene, const Ray &_ray) const {
  Spectrum tr(1.f);
  Ray ray(_ray);
  Point3f dest = _ray.at(_ray.tFar);
  do {
    auto itsOpt = scene.rayIntersect(ray);
    if (itsOpt.has_value()) {
      Intersection its = itsOpt.value();
      auto bsdf = its.shape->material->computeBSDF(its);
      if (!bsdf) {
        ray = Ray(its.position, dest);
        continue;
      } else {
        tr *= .0f;
      }
    }
    break;
  } while (1);
  return tr;
}

Spectrum VolPathIntegrator::li(const Ray &_ray, const Scene &scene,
                               std::shared_ptr<Sampler> sampler) const {
  Spectrum spectrum(.0f), beta(1.f);
  Ray ray(_ray);
  int step = 0;

  float pdfPrev;
  bool isDeltaPrev = true;

  do {
    auto itsOpt = scene.rayIntersect(ray);

    float tmax = itsOpt.has_value() ? itsOpt->distance : FLT_MAX;
    float t;

    MediumIntersection mits;
    if (auto medium = ray.medium; medium)
      beta *= medium->sample(ray, tmax, &mits);

    // Sample a medium intersection
    if (mits.valid) {
      const Medium *medium = mits.medium;
      // Sample infiniteLights
      for (auto light : scene.infiniteLights) {
        auto lightSampleResult = light->sample(mits, sampler->next2D());
        Ray shadowRay(mits.position, lightSampleResult.direction, 1e-4f,
                      lightSampleResult.distance);
        Spectrum tr = Tr(scene, shadowRay),
                 f = medium->phase->f(-ray.direction, shadowRay.direction);
        if (!tr.isZero() && !f.isZero()) {
          float pdf = convertPDF(lightSampleResult, mits);
          Spectrum emission = lightSampleResult.energy;
          spectrum += beta * f * tr * emission / pdf; // TODO MIS
        }
      }

      // Sample light in scene
      float pdfLight = .0f;
      auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
      if (light && pdfLight != .0f) {
        auto lightSampleResult = light->sample(mits, sampler->next2D());
        Ray shadowRay(mits.position, lightSampleResult.direction, 1e-4f,
                      lightSampleResult.distance);
        Spectrum tr = Tr(scene, shadowRay),
                 f = medium->phase->f(-ray.direction, shadowRay.direction);
        if (!tr.isZero() && !f.isZero()) {
          float pdf = convertPDF(lightSampleResult, mits);
          Spectrum emission = lightSampleResult.energy;
          spectrum += beta * f * tr * emission / pdf; // TODO MIS
        }
      }

      // Sample phase function to spwan the ray
      PhaseSampleResult phaseSampleResult =
          medium->phase->sample(-ray.direction, sampler->next2D());

      beta *= phaseSampleResult.weight;
      if (beta.isZero())
        break;
      ray = Ray(mits.position, phaseSampleResult.wi, 1e-4f, FLT_MAX);

      pdfPrev = phaseSampleResult.pdf;
      isDeltaPrev = false;
      ++step;
    }
    // Escape the scene
    else if (!itsOpt.has_value()) {
      for (auto light : scene.infiniteLights) {
        float pdf = light->pdf(ray);
        float misw = isDeltaPrev ? 1.f : powerHeuristic(pdfPrev, pdf);
        spectrum += beta * misw * light->evaluateEmission(ray);
      }
      break;
    }
    // Through the medium hit the surface
    else {
      Intersection its = itsOpt.value();
      auto bsdf = its.shape->material->computeBSDF(its);

      if (!bsdf) {
        ray = Ray(its.position, ray.direction, 1e-4f, FLT_MAX);
        continue;
      }

      // Hit the light in scene
      if (auto light = its.shape->light; light) {
        float pdf = light->pdf(ray, its);
        pdf *= scene.pdf(light);
        float misw = isDeltaPrev ? 1.f : powerHeuristic(pdfPrev, pdf);
        spectrum += beta * misw * light->evaluateEmission(its, -ray.direction);
      }

      if (step > 3)
        break;

      // Sample the Infinite light
      for (auto light : scene.infiniteLights) {
        auto lightSampleResult = light->sample(its, sampler->next2D());
        Ray shadowRay(its.position, lightSampleResult.direction, 1e-4f,
                      lightSampleResult.distance);
        Spectrum tr = Tr(scene, shadowRay),
                 f = bsdf->f(-ray.direction, shadowRay.direction);
        float pdf = convertPDF(lightSampleResult, its);
        if (!tr.isZero() && !f.isZero()) {
          float misw =
              lightSampleResult.isDelta
                  ? 1.f
                  : powerHeuristic(
                        pdf, bsdf->pdf(-ray.direction, shadowRay.direction));
          spectrum += beta * tr * misw * f * lightSampleResult.energy / pdf;
        }
      }

      // Sample the light in scene
      float pdfLight = .0f;
      auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
      if (light && pdfLight != .0f) {
        auto lightSampleResult = light->sample(its, sampler->next2D());
        Ray shadowRay(its.position, lightSampleResult.direction, 1e-4f,
                      lightSampleResult.distance);
        Spectrum tr = Tr(scene, shadowRay),
                 f = bsdf->f(-ray.direction, shadowRay.direction);
        lightSampleResult.pdf *= pdfLight;
        float pdf = convertPDF(lightSampleResult, its);
        if (!tr.isZero() && !f.isZero()) {
          float misw =
              lightSampleResult.isDelta
                  ? 1.f
                  : powerHeuristic(
                        pdf, bsdf->pdf(-ray.direction, shadowRay.direction));
          spectrum += beta * tr * misw * f * lightSampleResult.energy / pdf;
        }
      }

      // Sample the bsdf spwan the ray
      auto bsdfSampleResult = bsdf->sample(-ray.direction, sampler->next2D());
      beta *= bsdfSampleResult.weight;
      if (beta.isZero())
        break;
      ray = Ray(its.position, bsdfSampleResult.wi, 1e-4f, FLT_MAX);
      pdfPrev = bsdfSampleResult.pdf;
      isDeltaPrev = bsdfSampleResult.type == BSDFType::Specular;
      ++step;
    }

  } while (1);

  return spectrum;
}

REGISTER_CLASS(VolPathIntegrator, "volpath")