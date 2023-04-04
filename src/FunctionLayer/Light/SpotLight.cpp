#include "SpotLight.h"
#include <ResourceLayer/Factory.h>

SpotLight::SpotLight(const Json &json) : Light(json) {
  position = fetchRequired<Point3f>(json, "position");
  energy = fetchRequired<Spectrum>(json, "energy");
  type = LightType::SpotLight;
}

//! 由于点光源不会与光线发生相交，故该函数实际上不会被调用
Spectrum SpotLight::evaluateEmission(const Intersection &intersection,
                                     const Vector3f &wo) const {
  return Spectrum(.0f);
}

//! 做MIS的前提是其他采样方法采样到了该光源，而点光源只能由直接采样的方法得到，因此该函数也不会被调用
float SpotLight::pdf(const Ray &ray, const Intersection &intersection) const {
  return .0f;
}

LightSampleResult SpotLight::sample(const Intersection &shadingPoint,
                                    const Vector2f &sample) const {
  Vector3f shadingPoint2sample = position - shadingPoint.position;
  return LightSampleResult{energy,
                           normalize(shadingPoint2sample),
                           shadingPoint2sample.length() - EPSILON,
                           Vector3f(),
                           1.f,
                           true,
                           type};
}

void SpotLight::sampleLe(Vector2f u_position, Vector2f u_direction,
                         Ray *photonRay, float *pdf, Spectrum *Le,
                         Vector3f *lightNormal) const { // No implementation now
  // TODO
}

REGISTER_CLASS(SpotLight, "spotLight")