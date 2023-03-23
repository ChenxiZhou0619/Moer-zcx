#pragma once

#include "../Integrator.h"

class VolumetricPathTracer : public PixelIntegrator {
public:
  VolumetricPathTracer() = delete;

  VolumetricPathTracer(const Json &json);

  virtual ~VolumetricPathTracer() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override = 0;

protected:
  Spectrum Transmittance(
      const Scene &scene, Ray ray,
      std::function<Spectrum(const Medium *, Ray, float)> tr_tracker) const;

  std::optional<Intersection> TransmittanceRayIntersect(
      const Scene &scene, Ray ray, Spectrum *Tr,
      std::function<Spectrum(const Medium *, Ray, float)> tr_tracker) const;

  void setRayMedium(Vector3f direction, Vector3f normal,
                    std::shared_ptr<Material> material, Ray *ray) const;
};