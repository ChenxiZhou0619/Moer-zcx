#pragma once
#include "Material.h"

class EmptyMaterial : public Material {
public:
  EmptyMaterial() = default;

  EmptyMaterial(const Json &json) {
    // do nothing
  }

  // Return nullptr for empty material
  virtual std::shared_ptr<BSDF>
  computeBSDF(const Intersection &intersection) const override {
    return nullptr;
  }
};