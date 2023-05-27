#pragma once
#include "Material.h"

class RoughConductor : public Material {
public:
  RoughConductor();

  RoughConductor(const Json &json);

  virtual std::shared_ptr<BSDF>
  computeBSDF(const Intersection &intersection) const override;

private:
  std::shared_ptr<Texture<Spectrum>> albedo;
  float roughness;
};