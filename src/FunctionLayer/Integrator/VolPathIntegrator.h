#pragma once

#include "Integrator.h"

class VolPathIntegrator : public Integrator {
public:
  VolPathIntegrator() = default;

  VolPathIntegrator(const Json &json) : Integrator(json) {
    //
  }

  virtual ~VolPathIntegrator() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

protected:
  Spectrum Tr(const Scene &scene, const Ray &ray) const;

  void scatterSurface(Vector3f direction, Vector3f normal,
                      std::shared_ptr<Material> material, Ray *ray) const;
};