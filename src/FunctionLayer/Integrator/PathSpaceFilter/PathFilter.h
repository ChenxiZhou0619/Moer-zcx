#pragma once
#include "../Integrator.h"
#include <mutex>
struct ShadingVertex {
public:
  Point3f position;
  Vector3f normal;

  Spectrum weight;
  Spectrum Li = Spectrum(.0f), Li_refine = Spectrum(.0f);
  Spectrum Le = Spectrum(.0f);

  bool valid = false; /* Whether valid for filtering Li */
  Vector2i pixelLoc;  /* Pixel location on film */
  float pixelSA;      /* Pixel Solid angle to the focalpoint */
  float distSqr;
};

//* Store all shading vertices in scene
struct SVGroup {
public:
  SVGroup() = default;

  SVGroup(size_t size) { svs.reserve(size); }

  void emplace_back(ShadingVertex sv) {
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

  std::vector<ShadingVertex> svs;

  std::mutex mtx;
};

class PathFilter : public Integrator {
public:
  PathFilter() = default;

  virtual ~PathFilter() = default;

  PathFilter(const Json &json) : Integrator(json) {
    maxDepth = fetchOptional(json, "maxPathLength", 5);
    iteration = fetchOptional(json, "iteration", 3);
  }

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override {

    // No implementation
  }

  virtual void render(const Camera &camera, const Scene &scene,
                      std::shared_ptr<Sampler> sampler, int spp) const override;

  void samplePath(Ray ray, const Scene &scene, std::shared_ptr<Sampler>,
                  ShadingVertex *sv) const;

private:
  int maxDepth = 5;
  int iteration = 3;
};