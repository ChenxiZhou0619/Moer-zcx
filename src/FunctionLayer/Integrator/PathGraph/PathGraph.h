#pragma once

#include "../Integrator.h"
#include "Vertex.h"

class PathGraphIntegrator : public Integrator {
public:
  PathGraphIntegrator() = default;

  PathGraphIntegrator(const Json &json) : Integrator(json) {
    maxPathLength = fetchOptional(json, "maxPathLength", 5);
  }

  virtual ~PathGraphIntegrator() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

protected:
  int maxPathLength = 5;
  int rrThresholdLength = 5;
  float rrThreshold = .95f;
};
