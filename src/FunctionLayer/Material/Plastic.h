#pragma once
#include "Material.h"

class Plastic : public Material {
public:
  Plastic();

  Plastic(const Json &json);

  virtual std::shared_ptr<BSDF>
  computeBSDF(const Intersection &intersection) const override;

private:
  std::shared_ptr<Texture<Spectrum>> kd;
  float roughness, eta;
};