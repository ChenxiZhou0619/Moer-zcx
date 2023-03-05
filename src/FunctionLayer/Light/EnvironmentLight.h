#pragma once
#include "Light.h"
#include <CoreLayer/Math/Distribution.h>
#include <FunctionLayer/Ray/Ray.h>
#include <FunctionLayer/Texture/ImageTexture.h>
class EnvironmentLight : public InfiniteLight {
public:
  EnvironmentLight() = delete;

  EnvironmentLight(const Json &json);

  virtual Spectrum evaluateEmission(const Ray &ray) const override;

  virtual float pdf(const Ray &ray) const override;

  virtual float pdf(const Ray &ray,
                    const Intersection &its) const override final {
    // This will not be invoke
    return .0f;
  }

  virtual LightSampleResult sample(const Intersection &shadingPoint,
                                   const Vector2f &sample) const override;

private:
  std::shared_ptr<Texture<Spectrum>> environmentMap;
  Distribution1D<Vector2i> energyDistribution;
};