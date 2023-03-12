#pragma once

#include "../Integrator.h"
#include "Vertex.h"

class PathGraphIntegrator : public Integrator {
public:
  PathGraphIntegrator() = default;

  PathGraphIntegrator(const Json &json) : Integrator(json) {
    maxPathLength = fetchOptional(json, "maxPathLength", 5);
  }

  virtual ~PathGraphIntegrator() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

  LightPath constructPath(Vector2i pixelLoc, const Ray &ray, const Scene &scene,
                          std::shared_ptr<Sampler> sampler) const;

  virtual void render(const Camera &camera, const Scene &scene,
                      std::shared_ptr<Sampler> sampler, int spp) const override;

protected:
  int maxPathLength = 5;
  int rrThresholdLength = 5;
  float rrThreshold = .95f;
};
