#pragma once

#include "Integrator.h"

class PathIntegrator : public Integrator {
public:
  PathIntegrator() = default;

  PathIntegrator(const Json &json) : Integrator(json) {
    //
  }

  virtual ~PathIntegrator() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;
};