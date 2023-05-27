#include "Plastic.h"
#include "./BxDF/Microfacet.h"
#include <FunctionLayer/Texture/ConstantTexture.h>

Plastic::Plastic() {
  // default
  kd = std::make_shared<ConstantTexture<Spectrum>>(.5f);
  eta = 1.33f;
  roughness = 0.2;
}

Plastic::Plastic(const Json &json) : Material(json) {
  if (json["kd"].is_object()) {
    //* kd是一个json对象（表示纹理）
    kd = Factory::construct_class<Texture<Spectrum>>(json["kd"]);
  } else if (json["kd"].is_array()) {
    //* kd是一个json数组（表示常量kd）
    auto s = fetchRequired<Spectrum>(json, "kd");
    kd = std::make_shared<ConstantTexture<Spectrum>>(s);
  } else {
    std::cerr << "Error in kd format!\n";
    exit(1);
  }

  roughness = fetchOptional(json, "alpha", 0.2);
  eta = 1.33f;
}

std::shared_ptr<BSDF>
Plastic::computeBSDF(const Intersection &intersection) const {
  Vector3f normal, tangent, bitangent;
  computeShadingGeometry(intersection, &normal, &tangent, &bitangent);

  Spectrum s = kd->evaluate(intersection);
  return std::make_shared<CookTorranceReflection>(normal, tangent, bitangent,
                                                  roughness, eta, s);
}

REGISTER_CLASS(Plastic, "plastic")