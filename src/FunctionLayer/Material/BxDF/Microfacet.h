#pragma once
#include "BSDF.h"

class GGXDistribution {
public:
  static float D(Vector3f normal, float alpha) {
    float a = PI * alpha * alpha;
    float b = (normal[0] / alpha) * (normal[0] / alpha) +
              (normal[2] / alpha) * (normal[2] / alpha) + normal[1] * normal[1];
    return 1.f / (a * b * b);
  }

  static float G1(Vector3f w, float alpha) {
    auto lambda = [&](Vector3f direction) {
      float a =
          -1 + std::sqrt(1.f + (alpha * alpha * direction[0] * direction[0] +
                                alpha * alpha * direction[2] * direction[2]) /
                                   (direction[1] * direction[1]));
      return a * .5f;
    };

    float a = 1.f / (1.f + lambda(w));
    return a;
  }

  static Vector3f sampleNormal(Vector3f v, float alpha, Vector2f sample) {

    Vector3f wh = normalize(Vector3f{v[0] * alpha, v[1], v[2] * alpha});
    float lensq = v[0] * v[0] + v[2] * v[2];
    Vector3f T1 = lensq > .0f ? Vector3f(-v[2], 0, v[0]) / std::sqrt(lensq)
                              : Vector3f(1, 0, 0),
             T2 = cross(wh, T1);
    float r = std::sqrt(sample[0]), phi = 2 * PI * sample[1];
    float t1 = r * std::cos(phi), t2 = r * std::sin(phi);
    float s = .5f * (1.f + wh[1]);
    t2 = (1.0f - s) * std::sqrt(1.0f - t1 * t1) + s * t2;

    Vector3f Nh = t1 * T1 + t2 * T2 +
                  std::sqrt(std::max(0.f, 1.f - t1 * t1 - t2 * t2)) * wh;
    Vector3f Ne =
        normalize(Vector3f(alpha * Nh[0], std::max(0.f, Nh[1]), alpha * Nh[2]));
    return Ne;
  }
};

class MicrofacetReflection : public BSDF {
public:
  MicrofacetReflection(const Vector3f &_normal, const Vector3f &_tangent,
                       const Vector3f &_bitangent, float _alpha,
                       Spectrum _albedo)
      : BSDF(_normal, _tangent, _bitangent, BSDFType::Specular), alpha(_alpha),
        albedo(_albedo) {}

  virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const override {
    Vector3f woLocal = toLocal(wo), wiLocal = toLocal(wi),
             nLocal = normalize(woLocal + wiLocal);
    if (woLocal[1] < .0f || wiLocal[1] < .0f || nLocal[1] < .0f)
      return .0f;
    float G = nDistrib.G1(woLocal, alpha) * nDistrib.G1(wiLocal, alpha);
    return albedo * G * nDistrib.D(nLocal, alpha) * .25f / woLocal[1];
  }

  virtual BSDFSampleResult sample(const Vector3f &wo,
                                  const Vector2f &sample) const override {
    Vector3f woLocal = toLocal(wo);
    Vector3f nLocal = nDistrib.sampleNormal(woLocal, alpha, sample);
    Vector3f wiLocal = Reflect(woLocal, nLocal), wi = toWorld(wiLocal);
    float pdfWi = nDistrib.G1(woLocal, alpha) * nDistrib.D(nLocal, alpha) *
                  .25f / woLocal[1];
    return {albedo * nDistrib.G1(wiLocal, alpha), wi, pdfWi,
            BSDFType::Specular};
  }

  virtual float pdf(const Vector3f &wo, const Vector3f &wi) const override {
    Vector3f woLocal = toLocal(wo), wiLocal = toLocal(wi),
             nLocal = normalize(woLocal + wiLocal);
    return nDistrib.G1(woLocal, alpha) * nDistrib.D(nLocal, alpha) * .25f /
           woLocal[1];
  }

private:
  float alpha;
  Spectrum albedo;
  GGXDistribution nDistrib;
};