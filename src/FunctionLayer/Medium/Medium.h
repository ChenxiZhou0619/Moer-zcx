
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

  virtual bool Sample_MajorantTracking(Ray ray, float tmax, Vector2f sample,
                                       MediumIntersection *mits, Spectrum *Tr,
                                       float *pdf) const = 0;

  virtual bool Sample_WeightedMajorant(Ray ray, float tmax, Vector2f sample,
                                       MediumIntersection *mits, Spectrum *Tr,
                                       float *pdf) const {
    // No implementation
  }

  virtual Spectrum Transmittance_RatioTracking(Ray ray, float t) const = 0;

  virtual Spectrum Transmittance_ResidualRatioTracking(Ray ray, float t) const {
    // No implementation
  }

public:
  std::shared_ptr<Phase> phase;
};