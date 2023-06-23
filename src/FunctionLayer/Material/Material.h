#pragma once
#include "./BSSRDF/BSSRDF.h"
#include "./BxDF/BSDF.h"
#include <CoreLayer/Math/Math.h>
#include <FunctionLayer/Medium/Medium.h>
#include <FunctionLayer/Shape/Intersection.h>
#include <FunctionLayer/Texture/NormalTexture.h>
#include <ResourceLayer/Factory.h>
#include <ResourceLayer/JsonUtil.h>
class Material {
public:
  Material() {
    // donothing
  }

  Material(const Json &json) {
    if (json.contains("normalmap"))
      normalMap = std::make_shared<NormalTexture>(json["normalmap"]);
    if (json.contains("medium"))
      medium = Factory::construct_class<Medium>(json["medium"]);
    // TODO set bssrdf
  }

  virtual std::shared_ptr<BSDF>
  computeBSDF(const Intersection &intersection) const = 0;

  void computeShadingGeometry(const Intersection &intersection,
                              Vector3f *normal, Vector3f *tangent,
                              Vector3f *bitangent) const;

  std::shared_ptr<BSSRDF> computeBSSRDF(const Intersection &intersection) const;

  const Medium *getMedium() const;

  void setMedium(std::shared_ptr<Medium> medium);

protected:
  std::shared_ptr<NormalTexture> normalMap;
  std::shared_ptr<Medium> medium;
  std::shared_ptr<BSSRDF> bssrdf;
};