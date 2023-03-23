#include "RegularTracking.h"
#include <FunctionLayer/Material/Material.h>
#include <FunctionLayer/Medium/Medium.h>
RegularTracking::RegularTracking(const Json &json) : PixelIntegrator(json) {
  maxDepth = fetchOptional(json, "maxPathLength", 5);
}

Spectrum RegularTracking::li(const Ray &_ray, const Scene &scene,
                             std::shared_ptr<Sampler> sampler) const {
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
      mediumInteraction = medium->Sample_RegularTracking(
          ray, tmax, sampler->next2D(), &mi, &Tr, &pdf);
      beta *= Tr / pdf;
    }

    if (mediumInteraction) {
      //* Handle medium interaction
      if (++depth > maxDepth)
        break;

      auto phase = mi.mp.phase;
      Spectrum sigma_a = mi.mp.sigma_a, sigma_s = mi.mp.sigma_s;

      //* 1. Sample Le
      for (auto light : scene.infiniteLights) {
        LightSampleResult result = light->sample(mi, sampler->next2D());

        Ray shadowRay(mi.position, result.direction, 1e-4f, result.distance);
        shadowRay.medium = medium;
        Vector3f wi = shadowRay.direction;

        Spectrum Tr = Transmittance(scene, shadowRay), p = phase->f(wo, wi);
        float pdf = convertPDF(result, mi),
              misw = result.isDelta ? 1.f
                                    : powerHeuristic(pdf, phase->pdf(wo, wi));

        if (!Tr.isZero() && !p.isZero()) {
          L += beta * sigma_s * p * Tr * result.energy * misw / pdf;
        }
      }

      float pdfLight = .0f;
      if (auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
          light && pdfLight != .0f) {
        LightSampleResult result = light->sample(mi, sampler->next2D());

        Ray shadowRay(mi.position, result.direction, 1e-4f, result.distance);
        shadowRay.medium = medium;
        Vector3f wi = shadowRay.direction;

        Spectrum Tr = Transmittance(scene, shadowRay), p = phase->f(wo, wi);
        float pdf = convertPDF(result, mi) * pdfLight,
              misw = result.isDelta ? 1.f
                                    : powerHeuristic(pdf, phase->pdf(wo, wi));

        if (!Tr.isZero() && !p.isZero()) {
          L += beta * sigma_s * p * Tr * result.energy * misw / pdf;
        }
      }

      //* 2. Sample fp
      {
        PhaseSampleResult result = phase->sample(wo, sampler->next2D());

        Ray shadowRay{mi.position, result.wi, 1e-4f, FLT_MAX};
        shadowRay.medium = medium;
        Vector3f wi = shadowRay.direction;

        Spectrum Tr(1.f), p = phase->f(wo, wi);
        auto Le_si = TransmittanceRayIntersect(scene, shadowRay, &Tr);

        if (Le_si && Le_si->shape->light) {
          auto light = Le_si->shape->light;
          float pdfLe = scene.pdf(light) * light->pdf(shadowRay, *Le_si),
                misw = result.isDelta ? 1.f : powerHeuristic(result.pdf, pdfLe);
          Spectrum energy =
              light->evaluateEmission(*Le_si, -shadowRay.direction);
          L += beta * sigma_s * result.weight * Tr * energy * misw;
        } else if (!Le_si) {
          //* Hit environtment
          for (auto light : scene.infiniteLights) {
            float pdf = result.pdf,
                  misw = result.isDelta
                             ? 1.f
                             : powerHeuristic(pdf, light->pdf(shadowRay));
            Spectrum energy = light->evaluateEmission(shadowRay);
            L += beta * sigma_s * result.weight * Tr * energy * misw;
          }
        }
      }

      PhaseSampleResult result = phase->sample(wo, sampler->next2D());

      ray = Ray{mi.position, result.wi, 1e-4f, FLT_MAX};
      ray.medium = medium;
      beta *= result.weight * sigma_s;

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

          Spectrum Tr = Transmittance(scene, shadowRay), f = bsdf->f(wo, wi);
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

          Spectrum Tr = Transmittance(scene, shadowRay), f = bsdf->f(wo, wi);
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
          auto Le_si = TransmittanceRayIntersect(scene, shadowRay, &Tr);
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
        continue;
      }
    }
  }

  return L;
}

/*
// For debug
Spectrum RegularTracking::li(const Ray &_ray, const Scene &scene,
                             std::shared_ptr<Sampler> sampler) const {
  auto si = scene.rayIntersect(_ray);

  Ray ray{si->position, _ray.direction, 1e-4f, FLT_MAX};
  setRayMedium(ray.direction, si->normal, si->shape->material, &ray);

  const Medium *medium = ray.medium;
  MediumIntersection mits;
  Spectrum Tr;
  float pdf;

  si = scene.rayIntersect(ray);
  bool inMedium = medium->Sample_RegularTracking(ray, si->distance, {.9f, .0f},
                                                 &mits, &Tr, &pdf);
  if (inMedium) {
    std::cout << "Yes!\n";

    std::cout << "Position : ";
    mits.position.debugPrint();

    std::cout << "Tr : ";
    Tr.debugPrint();

    float dist = (mits.position - ray.origin).length();

    Spectrum rtTR = medium->Transmittance_RegularTracking(ray, dist);

    std::cout << "rtTR : ";
    rtTR.debugPrint();

    exit(1);
  }
}
*/

Spectrum RegularTracking::Transmittance(const Scene &scene, Ray ray) const {
  Spectrum Tr(1.f);

  while (true) {
    const Medium *medium = ray.medium;

    auto si = scene.rayIntersect(ray);
    float dt = si ? si->distance : ray.tFar;

    if (medium)
      Tr *= medium->Transmittance_RegularTracking(ray, dt);

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

std::optional<Intersection>
RegularTracking::TransmittanceRayIntersect(const Scene &scene, Ray ray,
                                           Spectrum *Tr) const {
  float distance = .0f;

  while (true) {
    const Medium *medium = ray.medium;
    auto si = scene.rayIntersect(ray);

    float dt = si ? si->distance : ray.tFar;

    if (medium)
      *Tr *= medium->Transmittance_RegularTracking(ray, dt);

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

void RegularTracking::setRayMedium(Vector3f direction, Vector3f normal,
                                   std::shared_ptr<Material> material,
                                   Ray *ray) const {
  bool towardsInner = dot(direction, normal) < .0f;
  ray->medium = towardsInner ? material->getMedium() : nullptr;
}

// Just for debug
// svoid RegularTracking::render(const Camera &camera, const Scene &scene,
//                             std::shared_ptr<Sampler> sampler, int spp) const
//                             {
//  Ray ray = camera.sampleRayDifferentials(CameraSample(), {.5f, .5f});
//  Spectrum s = li(ray, scene, sampler);
//}

REGISTER_CLASS(RegularTracking, "regularTracking")