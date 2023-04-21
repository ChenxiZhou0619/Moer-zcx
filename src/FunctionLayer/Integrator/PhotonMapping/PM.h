#pragma once

#include "../Integrator.h"

class PhotonMapping : public Integrator {
public:
  PhotonMapping() = default;

  virtual ~PhotonMapping() = default;

  PhotonMapping(const Json &json);

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override {
    // No implementation
  }

  virtual void render(const Camera &camera, const Scene &scene,
                      std::shared_ptr<Sampler> sampler, int spp) const override;

protected:
  int maxDepth;
  int k;
  int N;
};