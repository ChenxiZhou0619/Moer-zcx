#include "DeltaTracking.h"

DeltaTracking::DeltaTracking(const Json &json) : VolumetricPathTracer(json) {
  maxDepth = fetchOptional(json, "maxPathLength", 5);
}

Spectrum DeltaTracking::li(const Ray &ray, const Scene &scene,
                           std::shared_ptr<Sampler> sampler) const {
  //
}