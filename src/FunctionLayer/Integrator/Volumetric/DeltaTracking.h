#pragma once
#include "VolumetricPathTracer.h"

class DeltaTracking : public VolumetricPathTracer {
public:
  DeltaTracking() = delete;

  DeltaTracking(const Json &json);

  virtual ~DeltaTracking() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

private:
  int maxDepth;
};