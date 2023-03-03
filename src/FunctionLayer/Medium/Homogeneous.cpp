#include "Homogeneous.h"
#include "IsotropicPhase.h"
#include <FunctionLayer/Shape/Intersection.h>
HomogeneousMedium::HomogeneousMedium(const Json &json) : Medium(json) {
  sigmaA = fetchRequired<Spectrum>(json, "sigmaA");
  sigmaS = fetchRequired<Spectrum>(json, "sigmaS");
  sigmaT = sigmaA[0] + sigmaS[0];
  phase = std::make_shared<IsotropicPhase>(); // TODO
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
  return tr;
}

REGISTER_CLASS(HomogeneousMedium, "homogeneous")