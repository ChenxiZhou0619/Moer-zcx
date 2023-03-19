
#pragma once
#include "./Phase/Phase.h"
#include <CoreLayer/ColorSpace/Spectrum.h>
#include <CoreLayer/Math/Math.h>
#include <FunctionLayer/Ray/Ray.h>
#include <ResourceLayer/Factory.h>
#include <ResourceLayer/JsonUtil.h>

struct MediumIntersection;

class Medium {
public:
  Medium() = delete;

  Medium(const Json &json);

  virtual Spectrum sample(const Ray &ray, float tmax, Vector2f sample,
                          MediumIntersection *mits) const = 0;

  virtual Spectrum Tr(Point3f origin, Vector3f direction,
                      float distance) const = 0;

  virtual Spectrum SigmaS(Point3f position) const = 0;

  //* New Functions

  virtual bool Sample_RegularTracking(Ray ray, float tmax, Vector2f sample,
                                      MediumIntersection *mits, Spectrum *Tr,
                                      float *pdf) const = 0;

  virtual Spectrum Transmittance_RegularTracking(Ray ray, float t) const = 0;

public:
  std::shared_ptr<Phase> phase;
};