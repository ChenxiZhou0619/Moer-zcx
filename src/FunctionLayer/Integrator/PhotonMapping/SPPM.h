#pragma once

//* Stochastic Progressive Photon Mapping

#include "../Integrator.h"

class SPPM : public Integrator {
public:
  SPPM() = default;

  virtual ~SPPM() = default;

  SPPM(const Json &json);

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override {
    // No implementation
  }

  virtual void render(const Camera &camera, const Scene &scene,
                      std::shared_ptr<Sampler> sampler, int spp) const override;

protected:
  // TODO
  int maxDepth = 5;
};