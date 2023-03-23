#pragma once
#include "VolumetricPathTracer.h"
class RegularTracking : public VolumetricPathTracer {
public:
  RegularTracking() = delete;

  RegularTracking(const Json &json);

  virtual ~RegularTracking() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

private:
  int maxDepth;
};