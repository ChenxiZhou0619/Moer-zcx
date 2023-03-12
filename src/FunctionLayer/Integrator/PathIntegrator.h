#pragma once

#include "Integrator.h"

class PathIntegrator : public PixelIntegrator {
public:
  PathIntegrator() = default;

  PathIntegrator(const Json &json) : PixelIntegrator(json) {
    maxPathLength = fetchOptional(json, "maxPathLength", 5);
  }

  virtual ~PathIntegrator() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

protected:
  int maxPathLength = 5;
  int rrThresholdLength = 5;
  float rrThreshold = .95f;
};