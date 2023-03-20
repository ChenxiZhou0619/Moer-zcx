#include "Homogeneous.h"
#include "./Phase/HGPhase.h"
#include "./Phase/IsotropicPhase.h"
#include <FunctionLayer/Shape/Intersection.h>
HomogeneousMedium::HomogeneousMedium(const Json &json) : Medium(json) {
  sigmaA = fetchRequired<Spectrum>(json, "sigmaA");
  sigmaS = fetchRequired<Spectrum>(json, "sigmaS");
  sigmaT = sigmaA[0] + sigmaS[0];
}

bool HomogeneousMedium::Sample_RegularTracking(Ray ray, float tmax,
                                               Vector2f sample,
                                               MediumIntersection *mits,
                                               Spectrum *Tr, float *pdf) const {
  float dt = -std::log(1 - sample[0]) / sigmaT;
  if (dt > tmax) {
    *Tr = Transmittance_RegularTracking(ray, tmax);
    *pdf = (*Tr)[0];
    return false;
  } else {
    mits->position = ray.at(dt);
    mits->mp.phase = phase;
    mits->mp.sigma_a = sigmaA;
    mits->mp.sigma_s = sigmaS;
    *Tr = Transmittance_RegularTracking(ray, dt);
    *pdf = (*Tr)[0] * sigmaT;
    return true;
  }
}

Spectrum HomogeneousMedium::Transmittance_RegularTracking(Ray ray,
                                                          float t) const {
  return Spectrum(std::exp(-t * sigmaT));
}

REGISTER_CLASS(HomogeneousMedium, "homogeneous")