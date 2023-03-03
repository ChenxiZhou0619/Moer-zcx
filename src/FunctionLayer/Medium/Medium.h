
#pragma once
#include <CoreLayer/ColorSpace/Spectrum.h>
#include <CoreLayer/Math/Math.h>
#include <FunctionLayer/Ray/Ray.h>
#include <ResourceLayer/Factory.h>
#include <ResourceLayer/JsonUtil.h>

struct MediumIntersection;

struct PhaseSampleResult {
  Spectrum weight;
  Vector3f wi;
  float pdf;
};

class Phase {
public:
  virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const = 0;

  virtual float pdf(const Vector3f &wo, const Vector3f &wi) const = 0;

  virtual PhaseSampleResult sample(const Vector3f &wo,
                                   Vector2f sample) const = 0;
};

class Medium {
public:
  Medium() = delete;

  Medium(const Json &json){};

  virtual Spectrum sample(const Ray &ray, float tmax,
                          MediumIntersection *mits) const = 0;

public:
  std::shared_ptr<Phase> phase;
};