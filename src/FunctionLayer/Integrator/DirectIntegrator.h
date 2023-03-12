#pragma once
#include "Integrator.h"

//* 只计算场景中的直接光照，从光线与几何体的交点采样光源并连接
class DirectIntegratorSampleLight : public PixelIntegrator {
public:
  DirectIntegratorSampleLight() = default;

  DirectIntegratorSampleLight(const Json &json) {}

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;
};

//* 只计算场景中的直接光照，从光线与几何体的交点采样BSDF延伸光线，看其是否与光源相交
class DirectIntegratorSampleBSDF : public PixelIntegrator {
public:
  DirectIntegratorSampleBSDF() = default;

  DirectIntegratorSampleBSDF(const Json &json) {}

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;
};