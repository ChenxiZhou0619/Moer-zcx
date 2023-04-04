#include "AreaLight.h"
#include <FunctionLayer/Material/BxDF/Warp.h>
#include <ResourceLayer/Factory.h>
AreaLight::AreaLight(const Json &json) : Light(json) {
  type = LightType::AreaLight;
  shape = Factory::construct_class<Shape>(json["shape"]);
  energy = fetchRequired<Spectrum>(json, "energy");
}

Spectrum AreaLight::evaluateEmission(const Intersection &intersection,
                                     const Vector3f &wo) const {
  bool oneSide = dot(wo, intersection.normal) > 0;
  return oneSide ? energy : Spectrum(.0f);
}

float AreaLight::pdf(const Ray &ray, const Intersection &intersection) const {
  float pdfdA = 1.f / shape->area;
  // Convert pdfdA to pdfdw
  float pdfdw = pdfdA * intersection.distance * intersection.distance /
                std::abs(dot(intersection.normal, ray.direction));
  return pdfdw;
}

LightSampleResult AreaLight::sample(const Intersection &shadingPoint,
                                    const Vector2f &sample) const {
  Intersection sampleResult;
  float pdf;
  shape->uniformSampleOnSurface(sample, &sampleResult, &pdf);
  Vector3f shadingPoint2sample = sampleResult.position - shadingPoint.position;

  bool oneSide = dot(normalize(shadingPoint2sample), sampleResult.normal) < 0;

  return {(oneSide ? energy : Spectrum(.0f)),
          normalize(shadingPoint2sample),
          shadingPoint2sample.length() - EPSILON,
          sampleResult.normal,
          pdf,
          false,
          type};
}

void AreaLight::sampleLe(Vector2f u_position, Vector2f u_direction,
                         Ray *photonRay, float *pdf, Spectrum *Le,
                         Vector3f *lightNormal) const {
  Intersection sampleResult;
  float pdf_position;
  shape->uniformSampleOnSurface(u_position, &sampleResult, &pdf_position);
  Vector3f dir_local = squareToCosineHemisphere(u_direction);
  float pdf_direction = squareToCosineHemispherePdf(dir_local);
  Vector3f dir_world = dir_local[0] * sampleResult.tangent +
                       dir_local[1] * sampleResult.normal +
                       dir_local[2] * sampleResult.bitangent;
  *photonRay = Ray{sampleResult.position, dir_world, 1e-4f, FLT_MAX};
  *pdf = pdf_position * pdf_direction;
  *Le = energy;
  *lightNormal = sampleResult.normal;
}

REGISTER_CLASS(AreaLight, "areaLight")