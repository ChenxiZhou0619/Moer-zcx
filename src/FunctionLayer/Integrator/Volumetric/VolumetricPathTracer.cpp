#include "VolumetricPathTracer.h"
#include <FunctionLayer/Material/Material.h>

VolumetricPathTracer::VolumetricPathTracer(const Json &json)
    : PixelIntegrator(json) {
  // do nothing
}

Spectrum VolumetricPathTracer::Transmittance(
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

std::optional<Intersection> VolumetricPathTracer::TransmittanceRayIntersect(
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

void VolumetricPathTracer::setRayMedium(Vector3f direction, Vector3f normal,
                                        std::shared_ptr<Material> material,
                                        Ray *ray) const {
  bool towardsInner = dot(direction, normal) < .0f;
  ray->medium = towardsInner ? material->getMedium() : nullptr;
}