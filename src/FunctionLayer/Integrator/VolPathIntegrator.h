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

  Spectrum Tr(const Scene &scene, const Ray &ray) const;
};