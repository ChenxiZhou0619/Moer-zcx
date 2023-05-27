#include "RoughConductor.h"
#include "./BxDF/Microfacet.h"
#include <FunctionLayer/Texture/ConstantTexture.h>
RoughConductor::RoughConductor() {
  // default
  albedo = std::make_shared<ConstantTexture<Spectrum>>(Spectrum(1.f));
}

RoughConductor::RoughConductor(const Json &json) : Material(json) {
  if (json["albedo"].is_object()) {
    //* albedo是一个json对象（表示纹理）
    albedo = Factory::construct_class<Texture<Spectrum>>(json["albedo"]);
  } else if (json["albedo"].is_array()) {
    //* albedo是一个json数组（表示常量albedo）
    auto s = fetchRequired<Spectrum>(json, "albedo");
    albedo = std::make_shared<ConstantTexture<Spectrum>>(s);
  } else {
    std::cerr << "Error in albedo format!\n";
    exit(1);
  }

  roughness = fetchOptional(json, "alpha", 0.5);
}

std::shared_ptr<BSDF>
RoughConductor::computeBSDF(const Intersection &intersection) const {
  Vector3f normal, tangent, bitangent;
  computeShadingGeometry(intersection, &normal, &tangent, &bitangent);

  Spectrum s = albedo->evaluate(intersection);
  return std::make_shared<MicrofacetReflection>(normal, tangent, bitangent,
                                                roughness, s);
}

REGISTER_CLASS(RoughConductor, "roughConductor")