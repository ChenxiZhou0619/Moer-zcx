#include "VolPathIntegrator.h"
#include <FunctionLayer/Material/Material.h>
#include <FunctionLayer/Medium/Medium.h>
#include <FunctionLayer/Medium/Phase/Phase.h>

Spectrum VolPathIntegrator::li(const Ray &_ray, const Scene &scene,
                               std::shared_ptr<Sampler> sampler) const {
  Spectrum spectrum(.0f), beta(1.f);
  Ray ray(_ray);
  int pathLength = 0;

  std::optional<Intersection> itsOpt = scene.rayIntersect(ray);
  bool hitSurface = itsOpt.has_value();
  do {
    const Medium *medium = ray.medium;
    float tmax = hitSurface ? itsOpt->distance : FLT_MAX;
    MediumIntersection mits;
    if (medium)
      beta *= medium->sample(ray, tmax, sampler->next2D(), &mits);

    // If sample a medium intersection
    if (mits.valid) {
      if (++pathLength > maxPathLength)
        break;

      Spectrum sigmaS = mits.sigmaS;
      // Next event estimation
      ///  1. Sample infinite lights
      for (auto light : scene.infiniteLights) {
        auto result = light->sample(mits, sampler->next2D());
        Ray shadowRay(mits.position, result.direction, 1e-4f, result.distance);
        shadowRay.medium = medium;
        Spectrum tr = Tr(scene, shadowRay),
                 f = medium->phase->f(-ray.direction, shadowRay.direction);
        if (!tr.isZero() && !f.isZero()) {
          float pdf = convertPDF(result, mits);
          float misw = result.isDelta
                           ? 1.f
                           : powerHeuristic(
                                 pdf, medium->phase->pdf(-ray.direction,
                                                         shadowRay.direction));
          spectrum += beta * sigmaS * f * tr * result.energy * misw / pdf;
        }
      }

      ///  2. Sample light in scene
      {
        float pdfLight = .0f;
        auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
        if (light && pdfLight != .0f) {
          auto result = light->sample(mits, sampler->next2D());
          Ray shadowRay(mits.position, result.direction, 1e-4f,
                        result.distance);
          shadowRay.medium = medium;
          Spectrum tr = Tr(scene, shadowRay),
                   f = medium->phase->f(-ray.direction, shadowRay.direction);
          if (!tr.isZero() && !f.isZero()) {
            float pdf = convertPDF(result, mits);
            float misw =
                result.isDelta
                    ? 1.f
                    : powerHeuristic(pdf,
                                     medium->phase->pdf(-ray.direction,
                                                        shadowRay.direction));
            spectrum += beta * sigmaS * f * tr * result.energy * misw / pdf;
          }
        }
      }

      // Sample the phase function
      {
        PhaseSampleResult result =
            medium->phase->sample(-ray.direction, sampler->next2D());
        beta *= result.weight * sigmaS;
        if (beta.isZero())
          break;
        ray = Ray(mits.position, result.wi, 1e-4f, FLT_MAX);
        ray.medium = medium;

        // Get the intersection of phase sampled ray and calculate radiance
        Spectrum tr(1.f);
        auto [first, final] = rayIntersectTr(scene, ray, &tr);
        if (!final.has_value()) {
          for (auto light : scene.infiniteLights) {
            float pdf = light->pdf(ray);
            float misw = result.isDelta ? 1.f : powerHeuristic(result.pdf, pdf);
            spectrum += beta * misw * tr * light->evaluateEmission(ray);
          }
        } else if (auto light = final->shape->light; light) {
          float pdf = scene.pdf(light);
          pdf *= light->pdf(ray, final.value());
          float misw = result.isDelta ? 1.f : powerHeuristic(result.pdf, pdf);
          spectrum += beta * misw * tr *
                      light->evaluateEmission(final.value(), -ray.direction);
        }
        // March the ray, update the status
        itsOpt = first;
        hitSurface = itsOpt.has_value();
      }
    } else {
      // All radiance which directly hitted by ray will be handled at the
      // end of the loop, except the ray generate by camera
      if (pathLength == 0) {
        // Check if the ray hit the light in scene
        if (hitSurface) {
          if (auto light = itsOpt->shape->light; light)
            spectrum +=
                beta * light->evaluateEmission(itsOpt.value(), -ray.direction);
        } else {
          // Hit the environment, terminate the loop
          for (auto light : scene.infiniteLights)
            spectrum += beta * light->evaluateEmission(ray);
          break;
        }
      }

      if (++pathLength > maxPathLength || !hitSurface)
        break;

      Intersection its = itsOpt.value();
      auto bsdf = its.shape->material->computeBSDF(its);
      // Just through this surface if hit empty material
      if (!bsdf) {
        --pathLength;
        ray = Ray(its.position + 1e-3f * ray.direction, ray.direction, 1e-4f,
                  FLT_MAX);
        setRayMedium(ray.direction, its.normal, its.shape->material, &ray);
        itsOpt = scene.rayIntersect(ray);
        hitSurface = itsOpt.has_value();
        continue;
      }

      // Next event estimation
      /// 1. Sample infinite lights
      for (auto light : scene.infiniteLights) {
        auto result = light->sample(its, sampler->next2D());
        Ray shadowRay(its.position, result.direction, 1e-4f, result.distance);
        setRayMedium(shadowRay.direction, its.normal, its.shape->material,
                     &shadowRay);
        Spectrum tr = Tr(scene, shadowRay),
                 f = bsdf->f(-ray.direction, shadowRay.direction);
        if (!tr.isZero() && !f.isZero()) {
          float pdf = convertPDF(result, its);
          float misw =
              result.isDelta
                  ? 1.f
                  : powerHeuristic(result.pdf, bsdf->pdf(-ray.direction,
                                                         shadowRay.direction));
          spectrum += beta * f * tr * result.energy * misw / pdf;
        }
      }

      /// 2. Sample light in scene
      {
        float pdfLight = .0f;
        auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
        if (light && pdfLight != .0f) {
          auto result = light->sample(its, sampler->next2D());
          Ray shadowRay(its.position, result.direction, 1e-4f, result.distance);
          setRayMedium(shadowRay.direction, its.normal, its.shape->material,
                       &shadowRay);
          Spectrum tr = Tr(scene, shadowRay),
                   f = bsdf->f(-ray.direction, shadowRay.direction);
          if (!tr.isZero() && !f.isZero()) {
            result.pdf *= pdfLight;
            float pdf = convertPDF(result, its);
            float misw =
                result.isDelta
                    ? 1.f
                    : powerHeuristic(
                          pdf, bsdf->pdf(-ray.direction, shadowRay.direction));
            spectrum += beta * f * tr * result.energy * misw / pdf;
          }
        }
      }

      // Sample the bsdf
      {
        BSDFSampleResult result =
            bsdf->sample(-ray.direction, sampler->next2D());
        beta *= result.weight;
        if (beta.isZero())
          break;
        ray = Ray(its.position, result.wi, 1e-4f, FLT_MAX);
        setRayMedium(ray.direction, its.normal, its.shape->material, &ray);

        // Get the intersection of bsdf sampled ray and caculate the radiance
        Spectrum tr(1.f);
        auto [first, final] = rayIntersectTr(scene, ray, &tr);
        if (!final.has_value()) {
          for (auto light : scene.infiniteLights) {
            float pdf = light->pdf(ray);
            float misw = (result.type == BSDFType::Specular)
                             ? 1.f
                             : powerHeuristic(result.pdf, pdf);
            spectrum += beta * misw * light->evaluateEmission(ray);
          }
        } else if (auto light = final->shape->light; light) {
          float pdf = scene.pdf(light);
          pdf *= light->pdf(ray, final.value());
          float misw = (result.type == BSDFType::Specular)
                           ? 1.f
                           : powerHeuristic(result.pdf, pdf);
          spectrum += beta * misw *
                      light->evaluateEmission(final.value(), -ray.direction);
        }
        // March the ray, update the status
        itsOpt = first;
        hitSurface = itsOpt.has_value();
      }
    }

    if (pathLength > rrThresholdLength && beta.maxComponent() < rrThreshold) {
      float q = std::max(1.f - beta.maxComponent(), .05f);
      if (sampler->next1D() < q)
        break;
      beta /= 1 - q;
    }

  } while (1);
  return spectrum;
}

// TODO Robust it
Spectrum VolPathIntegrator::Tr(const Scene &scene, const Ray &_ray) const {
  Spectrum tr(1.f);
  Ray ray(_ray);
  Point3f dest = ray.at(ray.tFar);

  int step = 0;
  while (++step) {
    if (step % 200 == 0) {
      std::cout << "To much step in tr\n";
    }

    if (ray.tFar < 5e-3f)
      break;

    auto itsOpt = scene.rayIntersect(ray);
    bool hitSurface = itsOpt.has_value();

    // If block by surface
    if (hitSurface && itsOpt->shape->material->computeBSDF(itsOpt.value())) {
      tr = Spectrum(.0f);
      break;
    }

    if (!hitSurface && (ray.tFar > 1e5f && ray.medium))
      break;

    if (auto medium = ray.medium; medium) {
      float t = hitSurface ? itsOpt->distance : ray.tFar;
      tr *= medium->Tr(ray.origin, ray.direction, t);
    }

    if (!hitSurface)
      break;

    // spwan the ray
    ray = Ray(itsOpt->position, dest);
    setRayMedium(ray.direction, itsOpt->normal, itsOpt->shape->material, &ray);
  }

  return tr;
}

// Ignore all empty surface
// !Notice, ray.tFar always FLT_MAX
// !Notice, intersection.distance should be the total distance
std::pair<std::optional<Intersection>, std::optional<Intersection>>
VolPathIntegrator::rayIntersectTr(const Scene &scene, Ray ray,
                                  Spectrum *tr) const {
  std::optional<Intersection> first = std::nullopt, second = std::nullopt;
  float distance = .0f;
  *tr = Spectrum(1.f);
  do {
    second = scene.rayIntersect(ray);
    bool hitSurface = second.has_value();

    if (first == std::nullopt)
      first = second;

    if (!hitSurface)
      return {first, second};

    distance += second->distance;

    if (auto medium = ray.medium; medium) {
      *tr *= medium->Tr(ray.origin, ray.direction, second->distance);
    }

    if (second->shape->material->computeBSDF(second.value())) {
      second->distance = distance;
      break;
    }

    ray = Ray(second->position, ray.direction, 1e-4f, FLT_MAX);
    setRayMedium(ray.direction, second->normal, second->shape->material, &ray);
  } while (1);
  return {first, second};
}

void VolPathIntegrator::setRayMedium(Vector3f direction, Vector3f normal,
                                     std::shared_ptr<Material> material,
                                     Ray *ray) const {
  bool towardsInner = dot(direction, normal) < .0f;
  ray->medium = towardsInner ? material->getMedium() : nullptr;
}

REGISTER_CLASS(VolPathIntegrator, "volpath")