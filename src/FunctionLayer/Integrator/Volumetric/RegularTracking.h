#pragma once
#include "../Integrator.h"

class RegularTracking : public PixelIntegrator {
public:
  RegularTracking() = delete;

  RegularTracking(const Json &json);

  virtual ~RegularTracking() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override;

  // virtual void render(const Camera &camera, const Scene &scene,
  //                     std::shared_ptr<Sampler> sampler, int spp) const
  //                     override;

protected:
  Spectrum Transmittance(const Scene &scene, Ray ray) const;

  std::optional<Intersection>
  TransmittanceRayIntersect(const Scene &scene, Ray ray, Spectrum *Tr) const;

  void setRayMedium(Vector3f direction, Vector3f normal,
                    std::shared_ptr<Material> material, Ray *ray) const;

private:
  int maxDepth;
};