#pragma once

#include "Integrator.h"

class VolPathIntegrator : public PixelIntegrator {
public:
  VolPathIntegrator() = default;

  VolPathIntegrator(const Json &json) : PixelIntegrator(json) {
    maxPathLength = fetchOptional(json, "maxPathLength", 5);
  }

  virtual ~VolPathIntegrator() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

protected:
  Spectrum Tr(const Scene &scene, const Ray &ray) const;

  // First firstItsOpt
  // Second finalItsOpt
  std::pair<std::optional<Intersection>, std::optional<Intersection>>
  rayIntersectTr(const Scene &scene, Ray ray, Spectrum *tr) const;

  void setRayMedium(Vector3f direction, Vector3f normal,
                    std::shared_ptr<Material> material, Ray *ray) const;

protected:
  int maxPathLength = 5;
  int rrThresholdLength = 5;
  float rrThreshold = .95f;
};