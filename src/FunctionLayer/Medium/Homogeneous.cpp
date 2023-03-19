#include "Homogeneous.h"
#include "./Phase/HGPhase.h"
#include "./Phase/IsotropicPhase.h"
#include <FunctionLayer/Shape/Intersection.h>
HomogeneousMedium::HomogeneousMedium(const Json &json) : Medium(json) {
  sigmaA = fetchRequired<Spectrum>(json, "sigmaA");
  sigmaS = fetchRequired<Spectrum>(json, "sigmaS");
  sigmaT = sigmaA[0] + sigmaS[0];
}

Spectrum HomogeneousMedium::Tr(Point3f origin, Vector3f direction,
                               float distance) const {
  // Tr = exp(-sigmaT * distance)
  return fm::exp(-sigmaT * distance);
}

Spectrum HomogeneousMedium::SigmaS(Point3f position) const { return sigmaS; }

Spectrum HomogeneousMedium::sample(const Ray &ray, float tmax, Vector2f sample,
                                   MediumIntersection *mits) const {
  float t = -fm::log(1 - sample[0]) / sigmaT;
  Spectrum tr(1.f);

  // If sample a medium intersection
  if (t < tmax) {
    mits->medium = this;
    mits->valid = true;
    mits->position = ray.at(t);
    mits->distance = t;
    mits->sigmaS = SigmaS(mits->position);
    tr = Tr(ray.origin, ray.direction, t);
    for (int i = 0; i < 3; ++i) {
      mits->pdf += sigmaT * tr[i];
    }
    mits->pdf /= 3.f;
  }
  // Sample through the medium
  else {
    mits->medium = nullptr;
    mits->valid = false;
    mits->position = ray.at(tmax);
    mits->distance = tmax;
    tr = Tr(ray.origin, ray.direction, tmax);
    for (int i = 0; i < 3; ++i) {
      mits->pdf += tr[i];
    }
    mits->pdf /= 3.f;
  }
  return tr / mits->pdf;
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