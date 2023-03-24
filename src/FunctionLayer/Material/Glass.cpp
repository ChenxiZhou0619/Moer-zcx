#include "Glass.h"
#include "./BxDF/SpecularDielectric.h"
GlassMaterial::GlassMaterial() {
  // default
  eta = 1.33f;
}

GlassMaterial::GlassMaterial(const Json &json) : Material(json) {
  //
  eta = 1.33f;
}

std::shared_ptr<BSDF>
GlassMaterial::computeBSDF(const Intersection &intersection) const {
  Vector3f normal, tangent, bitangent;
  computeShadingGeometry(intersection, &normal, &tangent, &bitangent);
  return std::make_shared<SpecularDielectric>(normal, tangent, bitangent, eta);
}

REGISTER_CLASS(GlassMaterial, "glass")