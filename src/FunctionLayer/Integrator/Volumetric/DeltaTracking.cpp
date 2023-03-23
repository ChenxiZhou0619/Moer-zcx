#include "DeltaTracking.h"
#include <FunctionLayer/Material/Material.h>
#include <FunctionLayer/Medium/Medium.h>
#include <FunctionLayer/Shape/Intersection.h>

int SampleMode(float sigma_a, float sigma_s, float sigma_n, float uMode) {
  uMode *= (sigma_a + sigma_s + sigma_n);
  if (uMode < sigma_a)
    return 0; // Absorb
  else if (uMode < sigma_a + sigma_s)
    return 1; // Scatter
  else
    return 2; // Null-collision
}

DeltaTracking::DeltaTracking(const Json &json) : VolumetricPathTracer(json) {
  maxDepth = fetchOptional(json, "maxPathLength", 5);
}

Spectrum DeltaTracking::li(const Ray &_ray, const Scene &scene,
                           std::shared_ptr<Sampler> sampler) const {
  auto tr_tracker_regular = [](const Medium *medium, Ray ray, float tmax) {
    // return medium->Transmittance_RegularTracking(ray, tmax);
    return medium->Transmittance_RatioTracking(ray, tmax);
  };

  Spectrum L(.0f), beta(1.f);
  Ray ray(_ray);
  int depth = 0;

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
      beta *= Tr / pdf; //! Notice, Tr/pdf always 1
    }

    if (mediumInteraction) {
      //* Handle medium interaction
      if (++depth > maxDepth)
        break;

      auto phase = mi.mp.phase;
      // TODO handle chromatic medium
      float sigma_a = mi.mp.sigma_a[0], sigma_s = mi.mp.sigma_s[0],
            sigma_n = std::max(.0f, mi.mp.sigma_maj[0] - sigma_a - sigma_s);

      int mode = SampleMode(sigma_a, sigma_s, sigma_n, sampler->next1D());

      if (mode == 0 /* Absorb */) {
        // TODO Handle emission in medium
        // Just terminate now
        break;
      } else if (mode == 1 /* Scatter */) {
        //! Notice, sigma_s shouldn't be in L computing
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
            L += beta * p * Tr * result.energy * misw / pdf;
          }
        }

        float pdfLight = .0f;
        if (auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
            light && pdfLight != .0f) {
          LightSampleResult result = light->sample(mi, sampler->next2D());

          Ray shadowRay(mi.position, result.direction, 1e-4f, result.distance);
          shadowRay.medium = medium;
          Vector3f wi = shadowRay.direction;

          Spectrum Tr = Transmittance(scene, shadowRay, tr_tracker_regular),
                   p = phase->f(wo, wi);
          float pdf = convertPDF(result, mi) * pdfLight,
                misw = result.isDelta ? 1.f
                                      : powerHeuristic(pdf, phase->pdf(wo, wi));

          if (!Tr.isZero() && !p.isZero()) {
            L += beta * p * Tr * result.energy * misw / pdf;
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

          if (Le_si && Le_si->shape->light) {
            auto light = Le_si->shape->light;
            float pdfLe = scene.pdf(light) * light->pdf(shadowRay, *Le_si),
                  misw =
                      result.isDelta ? 1.f : powerHeuristic(result.pdf, pdfLe);
            Spectrum energy =
                light->evaluateEmission(*Le_si, -shadowRay.direction);
            L += beta * result.weight * Tr * energy * misw;
          } else if (!Le_si) {
            //* Hit environtment
            for (auto light : scene.infiniteLights) {
              float pdf = result.pdf,
                    misw = result.isDelta
                               ? 1.f
                               : powerHeuristic(pdf, light->pdf(shadowRay));
              Spectrum energy = light->evaluateEmission(shadowRay);
              L += beta * result.weight * Tr * energy * misw;
            }
          }
        }

        PhaseSampleResult result = phase->sample(wo, sampler->next2D());

        ray = Ray{mi.position, result.wi, 1e-4f, FLT_MAX};
        ray.medium = medium;
        beta *= result.weight;
      } else if (mode == 2 /* Null-collison */) {
        // Just continue the ray
        --depth;

        ray = Ray{mi.position, ray.direction, 1e-4f, FLT_MAX};
        ray.medium = medium;
      } else {
        std::cerr << "Shouldn't arrive here!\n";
        exit(1);
      }
    } else {
      //* Handle surface interaction
      if (depth == 0) {
        if (!si) {
          for (auto light : scene.infiniteLights)
            L += beta * light->evaluateEmission(ray);
        } else if (auto light = si->shape->light; light) {
          L += beta * light->evaluateEmission(*si, wo);
        }
      }

      if (!si)
        break;

      auto material = si->shape->material;
      auto bsdf = material->computeBSDF(*si);

      if (!bsdf) {
        ray = Ray(si->position, ray.direction, 1e-4f, FLT_MAX);
        setRayMedium(ray.direction, si->normal, material, &ray);
        continue;
      } else {
        //* Handle surface interaction
        if (++depth > maxDepth)
          break;

        //* 1. Sample Le
        for (auto light : scene.infiniteLights) {
          LightSampleResult result = light->sample(*si, sampler->next2D());

          Ray shadowRay{si->position, result.direction, 1e-4f, result.distance};
          Vector3f wi = shadowRay.direction;
          setRayMedium(wi, si->normal, material, &shadowRay);

          Spectrum Tr = Transmittance(scene, shadowRay, tr_tracker_regular),
                   f = bsdf->f(wo, wi);
          float pdf = convertPDF(result, *si),
                misw = result.isDelta ? 1.f
                                      : powerHeuristic(pdf, bsdf->pdf(wo, wi));
          if (!Tr.isZero() && !f.isZero()) {
            L += beta * f * Tr * result.energy * misw / pdf;
          }
        }

        float pdfLight = .0f;
        if (auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
            light && pdfLight != .0f) {
          LightSampleResult result = light->sample(*si, sampler->next2D());

          Ray shadowRay(si->position, result.direction, 1e-4f, result.distance);
          setRayMedium(shadowRay.direction, si->normal, material, &shadowRay);
          Vector3f wi = shadowRay.direction;

          Spectrum Tr = Transmittance(scene, shadowRay, tr_tracker_regular),
                   f = bsdf->f(wo, wi);
          float pdf = convertPDF(result, *si) * pdfLight,
                misw = result.isDelta ? 1.f
                                      : powerHeuristic(pdf, bsdf->pdf(wo, wi));

          if (!Tr.isZero() && !f.isZero()) {
            L += beta * f * Tr * result.energy * misw / pdf;
          }
        }

        //* 2. Sample bsdf
        {
          BSDFSampleResult result = bsdf->sample(wo, sampler->next2D());
          Ray shadowRay = Ray{si->position, result.wi, 1e-4f, FLT_MAX};
          Vector3f wi = shadowRay.direction;
          setRayMedium(wi, si->normal, material, &shadowRay);

          Spectrum Tr(1.f);
          auto Le_si = TransmittanceRayIntersect(scene, shadowRay, &Tr,
                                                 tr_tracker_regular);
          if (Le_si && Le_si->shape->light) {
            auto light = Le_si->shape->light;
            float pdfLe = scene.pdf(light) * light->pdf(shadowRay, *Le_si),
                  misw = (result.type == BSDFType::Specular)
                             ? 1.f
                             : powerHeuristic(result.pdf, pdfLe);
            Spectrum energy =
                light->evaluateEmission(*Le_si, -shadowRay.direction);
            L += beta * result.weight * Tr * energy * misw;

          } else if (!Le_si) {
            //* Hit environment
            for (auto light : scene.infiniteLights) {
              float pdf = result.pdf,
                    misw = (result.type == BSDFType::Specular)
                               ? 1.f
                               : powerHeuristic(pdf, light->pdf(shadowRay));
              Spectrum energy = light->evaluateEmission(shadowRay);
              L += beta * result.weight * Tr * energy * misw;
            }
          }
        }

        //* Spwan the ray
        BSDFSampleResult result = bsdf->sample(wo, sampler->next2D());
        beta *= result.weight;

        ray = Ray{si->position, result.wi, 1e-4f, FLT_MAX};
        setRayMedium(ray.direction, si->normal, material, &ray);
      }
    }
  }

  return L;
}

REGISTER_CLASS(DeltaTracking, "deltaTracking")