#pragma once

//* Volumetric Stochastic Progressive Photon Mapping

#include "../Integrator.h"

class VolSPPM : public Integrator {
public:
  VolSPPM() = default;

  virtual ~VolSPPM() = default;

  VolSPPM(const Json &json);

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override {
    // No implementation
  }

  virtual void render(const Camera &camera, const Scene &scene,
                      std::shared_ptr<Sampler> sampler, int spp) const override;

protected:
  void setRayMedium(Vector3f direction, Vector3f normal,
                    std::shared_ptr<Material> material, Ray *ray) const;

protected:
  // TODO
  int maxDepth;
  int photonsPerIteration;
  float searchRadius;
};