#pragma once
#include "VolumetricPathTracer.h"

class WeightedDeltaTracking : public VolumetricPathTracer {
public:
  WeightedDeltaTracking() = delete;

  WeightedDeltaTracking(const Json &json);

  virtual ~WeightedDeltaTracking() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

private:
  int maxDepth;
};