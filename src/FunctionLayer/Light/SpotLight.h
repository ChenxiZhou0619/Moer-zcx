#pragma once
#include "Light.h"

class SpotLight : public Light {
public:
  SpotLight() = delete;

  SpotLight(const Json &json);

  virtual Spectrum evaluateEmission(const Intersection &intersection,
                                    const Vector3f &wo) const override;

  virtual LightSampleResult sample(const Intersection &shadingPoint,
                                   const Vector2f &sample) const override;

  virtual float pdf(const Ray &ray,
                    const Intersection &intersection) const override;

  virtual void sampleLe(Vector2f u_position, Vector2f u_direction,
                        Ray *photonRay, float *pdf, Spectrum *Le,
                        Vector3f *lightNormal) const override;

private:
  Point3f position;
  Spectrum energy;
};