#pragma once
#include "../Integrator.h"

//* 法线可视化
class NormalIntegrator : public PixelIntegrator {
public:
  NormalIntegrator() = default;

  NormalIntegrator(const Json &json) : PixelIntegrator(json) {}

  virtual ~NormalIntegrator() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;
};
