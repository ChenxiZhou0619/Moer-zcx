#pragma once
#include "BSDF.h"

class SpecularDielectric : public BSDF {
public:
  SpecularDielectric(const Vector3f &_normal, const Vector3f &_tangent,
                     const Vector3f &_bitangent, float eta)
      : BSDF(_normal, _tangent, _bitangent), eta(eta) {}

  virtual BSDFSampleResult sample(const Vector3f &wo,
                                  const Vector2f &sample) const override {
    Vector3f woLocal = toLocal(wo);

    auto FresnelTerm = [](float eta_wi, float eta_wo, float cos_wi,
                          float cos_wo) {
      float R1 = (eta_wo * cos_wi - eta_wi * cos_wo) /
                 (eta_wo * cos_wi + eta_wi * cos_wo),
            R2 = (eta_wi * cos_wi - eta_wo * cos_wo) /
                 (eta_wi * cos_wi + eta_wo * cos_wo);
      return (R1 * R1 + R2 * R2) * .5f;
    };

    float eta_wo = eta, eta_wi = 1.f;
    float cos_wo = woLocal[1];
    if (cos_wo > .0f)
      std::swap(eta_wo, eta_wi);

    //* sin_wi * eta_wi = sin_wo * eta_wo
    float sin_wo = fm::sqrt(std::max(.0f, 1 - cos_wo * cos_wo)),
          sin_wi = sin_wo * eta_wo / eta_wi;

    if (sin_wi >= 1.f) {
      //* Totally internal reflection
      Vector3f wiLocal{-woLocal[0], woLocal[1], -woLocal[2]};
      return {Spectrum(1.f), toWorld(wiLocal), 1.f, BSDFType::Specular};
    } else {
      //* Sample reflection / refraction according to fresnel term
      float cos_wi = fm::sqrt(std::max(.0f, 1 - sin_wi * sin_wi));

      //* Compute fresnel
      float F = FresnelTerm(eta_wi, eta_wo, cos_wi, std::abs(cos_wo));

      if (sample[0] < F) {
        //* 1. Specular reflection
        Vector3f wiLocal{-woLocal[0], woLocal[1], -woLocal[2]};
        return {Spectrum(1.f), toWorld(wiLocal), 1.f, BSDFType::Specular};
      } else {
        //* 2. Specular refraction
        float scale = eta_wo / eta_wi;
        if (cos_wi * cos_wo > 0)
          cos_wi *= -1;
        Vector3f wiLocal{-woLocal[0] * scale, cos_wi, -woLocal[2] * scale};
        float energyScale = (eta_wo / eta_wi) * (eta_wo / eta_wi);
        return {Spectrum(energyScale), toWorld(wiLocal), 1.f,
                BSDFType::Specular};
      }
    }
  }

  virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const override {
    return Spectrum(.0f);
  }

  virtual float pdf(const Vector3f &wo, const Vector3f &wi) const override {
    return .0f;
  }

private:
  float eta;
};