#include "Homogeneous.h"
#include "../Phase/HGPhase.h"
#include "../Phase/IsotropicPhase.h"
#include <FunctionLayer/Shape/Intersection.h>
HomogeneousMedium::HomogeneousMedium(const Json &json) : Medium(json) {
  sigmaA = fetchRequired<Spectrum>(json, "sigmaA");
  sigmaS = fetchRequired<Spectrum>(json, "sigmaS");
  sigmaT = sigmaA + sigmaS;
}

bool HomogeneousMedium::Sample_RegularTracking(Ray ray, float tmax,
                                               Vector2f sample,
                                               MediumIntersection *mits,
                                               Spectrum *Tr, float *pdf) const {

  auto sampleChannel = [&](float u) { return clamp<int>((int)(u * 3), 0, 2); };
  // random pick one channel
  int channel = sampleChannel(sample[1]);

  float dt = -std::log(1 - sample[0]) / sigmaT[channel];

  if (dt > tmax) {
    *Tr = Transmittance_RegularTracking(ray, tmax);
    *pdf = .0f;
    *pdf += ((*Tr)[0] + (*Tr)[1] + (*Tr)[2]) / 3.f;
    return false;
  } else {
    mits->position = ray.at(dt);
    mits->mp.phase = phase;
    mits->mp.sigma_a = sigmaA;
    mits->mp.sigma_s = sigmaS;
    mits->mp.sigma_maj = sigmaA + sigmaS;
    *Tr = Transmittance_RegularTracking(ray, dt);
    *pdf = .0f;
    *pdf +=
        ((*Tr)[0] * sigmaT[0] + (*Tr)[1] * sigmaT[1] + (*Tr)[2] * sigmaT[2]) /
        3.f;
    return true;
  }
}

Spectrum HomogeneousMedium::Transmittance_RegularTracking(Ray ray,
                                                          float t) const {
  return Spectrum(std::exp(-t * sigmaT[0]), std::exp(-t * sigmaT[1]),
                  std::exp(-t * sigmaT[2]));
}

//! Notice, Tr / pdf == 1
bool HomogeneousMedium::Sample_MajorantTracking(Ray ray, float tmax,
                                                Vector2f sample,
                                                MediumIntersection *mits,
                                                Spectrum *Tr,
                                                float *pdf) const {
  auto sampleChannel = [&](float u) { return clamp<int>((int)(u * 3), 0, 2); };
  // random pick one channel
  int channel = sampleChannel(sample[1]);
  float dt = -std::log(1 - sample[0]) / sigmaT[channel];
  if (dt > tmax) {
    *Tr = Spectrum(1.f);
    *pdf = 1.f;
    return false;
  } else {
    mits->position = ray.at(dt);
    mits->mp.phase = phase;
    mits->mp.sigma_a = sigmaA;
    mits->mp.sigma_s = sigmaS;
    *Tr = Spectrum(1.f);
    *pdf = 1.f;
    return true;
  }
}

Spectrum HomogeneousMedium::Transmittance_RatioTracking(Ray ray,
                                                        float t) const {
  return Spectrum(std::exp(-t * sigmaT[0]), std::exp(-t * sigmaT[1]),
                  std::exp(-t * sigmaT[2]));
}

REGISTER_CLASS(HomogeneousMedium, "homogeneous")