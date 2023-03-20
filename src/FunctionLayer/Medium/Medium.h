
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

  virtual bool Sample_RegularTracking(Ray ray, float tmax, Vector2f sample,
                                      MediumIntersection *mits, Spectrum *Tr,
                                      float *pdf) const = 0;

  virtual Spectrum Transmittance_RegularTracking(Ray ray, float t) const = 0;

public:
  std::shared_ptr<Phase> phase;
};