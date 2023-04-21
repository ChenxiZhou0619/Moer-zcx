#include "PM.h"

PhotonMapping::PhotonMapping(const Json &json) {
  maxDepth = fetchRequired<int>(json, "maxPathLength");
  k = fetchRequired<int>(json, "k");
  N = fetchRequired<int>(json, "N");
}

void PhotonMapping::render(const Camera &camera, const Scene &scene,
                           std::shared_ptr<Sampler> sampler, int spp) const {
  //* Generate N photons and store them in kd-tree
}