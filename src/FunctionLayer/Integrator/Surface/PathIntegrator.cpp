// PathIntegrator with mis

#include "PathIntegrator.h"
#include <FunctionLayer/Material/Material.h>

Spectrum PathIntegrator::li(const Ray &_ray, const Scene &scene,
                            std::shared_ptr<Sampler> sampler) const {
  Spectrum Li(.0f), beta(1.f);
  Ray ray(_ray);
  int bounces = 0;

  while (true) {
    auto itsOpt = scene.rayIntersect(ray);

    if (bounces == 0) {
      if (!itsOpt)
        for (auto light : scene.infiniteLights)
          Li += beta * light->evaluateEmission(ray);
      else if (auto light = itsOpt->shape->light; light)
        Li += beta * light->evaluateEmission(*itsOpt, -ray.direction);
    }

    if (!itsOpt || maxPathLength == 0)
      break;

    auto bsdf = itsOpt->shape->material->computeBSDF(*itsOpt);

    // Sample direct light
    Li += beta * sampleDirect(scene, *itsOpt, -ray.direction, bsdf, sampler);

    // Sample a direction according to bsdf to spwan the ray
    auto bsdfSampleResult = bsdf->sample(-ray.direction, sampler->next2D());
    beta *= bsdfSampleResult.weight;

    bool isTransmission = dot(bsdfSampleResult.wi, itsOpt->normal) < .0f;
    // TODO consider the bssrdf
    if (beta.isZero() || ++bounces >= maxPathLength)
      break;

    ray = Ray{itsOpt->position, bsdfSampleResult.wi, 1e-4f, FLT_MAX};
  }
  return Li;
}

REGISTER_CLASS(PathIntegrator, "path")