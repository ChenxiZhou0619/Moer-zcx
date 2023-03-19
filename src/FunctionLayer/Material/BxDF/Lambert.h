#pragma once
#include "BSDF.h"
#include "Warp.h"
class LambertReflection : public BSDF {
public:
  LambertReflection(const Vector3f &_normal, const Vector3f &_tangent,
                    const Vector3f &_bitangent, Spectrum _albedo)
      : BSDF(_normal, _tangent, _bitangent), albedo(_albedo) {}

  virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const override {
    Vector3f woLocal = toLocal(wo), wiLocal = toLocal(wi);
    if (woLocal[1] <= .0f || wiLocal[1] <= .0f)
      return Spectrum(0.f);
    return albedo * INV_PI * wiLocal[1];
  }

  virtual BSDFSampleResult sample(const Vector3f &wo,
                                  const Vector2f &sample) const override {
    Vector3f wi = squareToCosineHemisphere(sample);
    float pdf = squareToCosineHemispherePdf(wi);
    Vector3f w = toWorld(wi);
    if (std::isnan(w[0]) || std::isnan(w[1]) || std::isnan(w[2])) {
      std::cout << "Stop!\n";
    }

    return {albedo, toWorld(wi), pdf, BSDFType::Diffuse};
  }

  virtual float pdf(const Vector3f &wo, const Vector3f &wi) const override {
    return squareToCosineHemispherePdf(toLocal(wi));
  }

private:
  Spectrum albedo;
};