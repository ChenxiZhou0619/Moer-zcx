#pragma once

#include "../Integrator.h"

class WhittedIntegrator : public PixelIntegrator {
public:
  WhittedIntegrator() = default;

  WhittedIntegrator(const Json &json) : PixelIntegrator(json) {}

  virtual ~WhittedIntegrator() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;
};