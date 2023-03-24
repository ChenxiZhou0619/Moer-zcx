#pragma once
#include "Material.h"

class GlassMaterial : public Material {
public:
  GlassMaterial();

  GlassMaterial(const Json &json);

  virtual std::shared_ptr<BSDF>
  computeBSDF(const Intersection &intersection) const override;

private:
  float eta;
};