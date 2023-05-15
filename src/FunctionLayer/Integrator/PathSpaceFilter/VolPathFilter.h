#pragma once

#include "../Integrator.h"
#include <mutex>

struct VolShadingVertex {
public:
  Point3f position;

  Spectrum weight;
  Spectrum Li = Spectrum(.0f), Li_refine = Spectrum(.0f);
  Spectrum Le = Spectrum(.0f);

  float sigma_t;

  bool valid = false; /* Whether valid for filtering Li */
  Vector2i pixelLoc;  /* Pixel location on film */
  float pixelSA;      /* Pixel Solid angle to the focalpoint */
  float distSqr;
};

//* Store all shading vertices in scene
struct VolSVGroup {
public:
  VolSVGroup() = default;

  VolSVGroup(size_t size) { svs.reserve(size); }

  void emplace_back(VolShadingVertex sv) {
    std::lock_guard<std::mutex> lk_gd(mtx);
    svs.emplace_back(sv);
  }

  inline size_t kdtree_get_point_count() const { return svs.size(); }

  inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
    return svs[idx].position[dim];
  }

  // Just return false
  template <class BBox> bool kdtree_get_bbox(BBox & /*bb*/) const {
    return false;
  }

  std::vector<VolShadingVertex> svs;

  std::mutex mtx;
};

class VolPathFilter : public Integrator {
public:
  VolPathFilter() = default;

  virtual ~VolPathFilter() = default;

  VolPathFilter(const Json &json) : Integrator(json) {
    maxDepth = fetchOptional(json, "maxPathLength", 5);
    iteration = fetchOptional(json, "iteration", 3);
    k = fetchRequired<int>(json, "k");
  }

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override {
    // No implementation
  }

  virtual void render(const Camera &camera, const Scene &scene,
                      std::shared_ptr<Sampler> sampler, int spp) const override;

  void samplePath(Ray ray, const Scene &scene, std::shared_ptr<Sampler> sampler,
                  VolShadingVertex *sv) const;

private:
  void setRayMedium(Vector3f direction, Vector3f normal,
                    std::shared_ptr<Material> material, Ray *ray) const;

  Spectrum Transmittance(
      const Scene &scene, Ray ray,
      std::function<Spectrum(const Medium *, Ray, float)> tr_tracker) const;

  std::optional<Intersection> TransmittanceRayIntersect(
      const Scene &scene, Ray ray, Spectrum *Tr,
      std::function<Spectrum(const Medium *, Ray, float)> tr_tracker) const;

  int maxDepth = 5;
  int iteration = 3;
  int k;
};