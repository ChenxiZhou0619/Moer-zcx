#include "Material.h"

std::shared_ptr<BSSRDF>
Material::computeBSSRDF(const Intersection &intersection) const {
  return bssrdf;
}

const Medium *Material::getMedium() const {
  if (medium)
    return medium.get();
  return nullptr;
}

void Material::setMedium(std::shared_ptr<Medium> _medium) { medium = _medium; }

void Material::computeShadingGeometry(const Intersection &intersection,
                                      Vector3f *normal, Vector3f *tangent,
                                      Vector3f *bitangent) const {
  if (!normalMap) {
    *normal = intersection.normal;
    *tangent = intersection.tangent;
    *bitangent = intersection.bitangent;
  } else {
    Vector3f localNormal = normalMap->evaluate(intersection);
    *normal = normalize(localNormal[0] * intersection.tangent +
                        localNormal[1] * intersection.bitangent +
                        localNormal[2] * intersection.normal);
    *tangent = intersection.tangent;
    *bitangent = normalize(cross(*tangent, *normal));
  }
}