#pragma once
#include "Phase.h"

class HenyeyGrennsteinPhase : public Phase {
public:
  HenyeyGrennsteinPhase(float _g) : g(_g) {}

  virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const override;

  virtual float pdf(const Vector3f &wo, const Vector3f &wi) const override;

  virtual PhaseSampleResult sample(const Vector3f &wo,
                                   Vector2f sample) const override;

public:
  float g = 0.5;
};
