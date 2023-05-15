/**
 * @file VolPathGraph.h
 * @author Chenxi Zhou
 * @brief This only support volume only scenes
 * @version 0.1
 * @date 2023-05-15
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
#include "../Integrator.h"
#include "Vertex.h"

struct VolumeShadingPoint {
  //* Scatter point properties
  Point3f position;             // position of shading point
  std::shared_ptr<Phase> phase; // phase function of the shading point

  //* Direct lighting properties
  Vector3f wo; // direction towards camera
  struct DirectLight {
    Vector3f wi;        // direction of direct light
    Spectrum li{.0f};   // incident light on wi
    float weight = .0f; // weight of the edge
  } nee, phs;           // direct light sample by two distributions

  //* Indirect lighting properties
  Vector3f wj;
  VolumeShadingPoint *xj = nullptr;
  Spectrum lj{.0f}, weight_j{.0f};
};

struct VSPGroup {
public:
  VSPGroup() = default;

  VSPGroup(size_t size) { points.reserve(size); }

  size_t emplace_back(VolumeShadingPoint p) {
    std::lock_guard<std::mutex> lkgd(mtx);
    points.emplace_back(p);
    return points.size() - 1;
  }

  inline size_t kdtree_get_point_count() { return points.size(); }

  inline float kdtree_get_pt(const size_t idx, const size_t dim) const {
    return points[idx].position[dim];
  }

  // Just return false
  template <class BBox> bool kdtree_get_bbox(BBox & /*bb*/) const {
    return false;
  }

  std::vector<VolumeShadingPoint> points;

  std::mutex mtx;
};

class VolPathGraph : public Integrator {
public:
  VolPathGraph() = default;

  VolPathGraph(const Json &json) : Integrator(json) {
    maxPathLength = fetchOptional(json, "maxPathLength", 5);
  }

  virtual ~VolPathGraph() = default;

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override {
    // No implementation
  }

  virtual void render(const Camera &camera, const Scene &scene,
                      std::shared_ptr<Sampler> sampler, int spp) const override;

  Spectrum constructPath(Ray ray, const Scene &scene,
                         std::shared_ptr<Sampler> sampler, VSPGroup *group,
                         int *firstShadingPointIdx) const;

private:
  void setRayMedium(Vector3f direction, Vector3f normal,
                    std::shared_ptr<Material> material, Ray *ray) const;

  Spectrum Transmittance(
      const Scene &scene, Ray ray,
      std::function<Spectrum(const Medium *, Ray, float)> tr_tracker) const;

  std::optional<Intersection> TransmittanceRayIntersect(
      const Scene &scene, Ray ray, Spectrum *Tr,
      std::function<Spectrum(const Medium *, Ray, float)> tr_tracker) const;

  int maxPathLength;
  int iterations = 1;
};