#pragma once
#include "Phase.h"
#include <FunctionLayer/Material/BxDF/Warp.h>

class IsotropicPhase : public Phase {
public:
  IsotropicPhase() = default;

  virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const override {
    return 0.25f * INV_PI;
  }

  virtual float pdf(const Vector3f &wo, const Vector3f &wi) const override {
    return 0.25f * INV_PI;
  }

  virtual PhaseSampleResult sample(const Vector3f &wo,
                                   Vector2f sample) const override {
    return PhaseSampleResult{1, squareToUniformSphere(sample), 0.25f * INV_PI};
  }
};