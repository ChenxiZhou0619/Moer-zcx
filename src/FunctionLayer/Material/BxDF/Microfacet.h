#pragma once
#include "./MicrofacetDistribution.h"
#include "BSDF.h"
#include "Warp.h"
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
  static GGXDistribution nDistrib;
};

class CookTorranceReflection : public BSDF {
public:
  CookTorranceReflection(const Vector3f &_normal, const Vector3f &_tangent,
                         const Vector3f &_bitangent, float _alpha, float _eta,
                         Spectrum _kd)
      : BSDF(_normal, _tangent, _bitangent, BSDFType::Diffuse), alpha(_alpha),
        eta(_eta), kd(_kd), ks(Spectrum(1.f) - _kd) {}

  virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const override {
    Vector3f woLocal = toLocal(wo), wiLocal = toLocal(wi),
             nLocal = normalize(woLocal + wiLocal);
    if (woLocal[1] < .0f || wiLocal[1] < .0f || nLocal[1] < .0f)
      return .0f;
    float G = nDistrib.G1(woLocal, alpha) * nDistrib.G1(wiLocal, alpha),
          F = FresnelSchlick(wiLocal, eta);
    Spectrum specularTerm =
        F * G * nDistrib.D(nLocal, alpha) * .25f / woLocal[1];

    Spectrum diffuseTerm = INV_PI * wiLocal[1];

    return kd * diffuseTerm + ks * specularTerm;
  }

  virtual BSDFSampleResult sample(const Vector3f &wo,
                                  const Vector2f &sample) const override {

    float diffuseProb = (kd[0] + kd[1] + kd[2]) / 3.f,
          specularProb = (ks[0] + ks[1] + ks[2]) / 3.f;

    Vector3f wi;
    if (sample[0] + sample[1] < diffuseProb) {
      // sample the diffuse
      Vector3f wiLocal = squareToCosineHemisphere(sample);
      wi = toWorld(wiLocal);
    } else {
      Vector3f woLocal = toLocal(wo);
      Vector3f nLocal = nDistrib.sampleNormal(woLocal, alpha, sample);
      Vector3f wiLocal = Reflect(woLocal, nLocal);
      wi = toWorld(wiLocal);
    }

    float pdf = this->pdf(wo, wi);
    Spectrum weight = f(wo, wi) / pdf;
    return {weight, wi, pdf, BSDFType::Diffuse};
  }

  virtual float pdf(const Vector3f &wo, const Vector3f &wi) const override {
    float diffuseProb = (kd[0] + kd[1] + kd[2]) / 3.f,
          specularProb = (ks[0] + ks[1] + ks[2]) / 3.f;
    Vector3f woLocal = toLocal(wo), wiLocal = toLocal(wi),
             nLocal = normalize(woLocal + wiLocal);
    float pdfSpecular = nDistrib.G1(woLocal, alpha) *
                        nDistrib.D(nLocal, alpha) * .25f / woLocal[1];
    float pdfDiffuse = squareToCosineHemispherePdf(wiLocal);
    return diffuseProb * pdfDiffuse + specularProb * pdfSpecular;
  }

private:
  float alpha, eta;
  Spectrum kd, ks;
  static GGXDistribution nDistrib;
};