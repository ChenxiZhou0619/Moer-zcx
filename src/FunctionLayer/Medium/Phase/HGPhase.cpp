#include "HGPhase.h"

Spectrum HenyeyGrennsteinPhase::f(const Vector3f &wo,
                                  const Vector3f &wi) const {
  float cosTheta = dot(wo, wi);
  float denom = 1 + g * g + 2 * g * cosTheta;
  return 0.25 * INV_PI * (1 - g * g) / (denom * fm::sqrt(denom));
}

float HenyeyGrennsteinPhase::pdf(const Vector3f &wo, const Vector3f &wi) const {
  float cosTheta = dot(wo, wi);
  float denom = 1 + g * g + 2 * g * cosTheta;
  return 0.25 * INV_PI * (1 - g * g) / (denom * fm::sqrt(denom));
}

PhaseSampleResult HenyeyGrennsteinPhase::sample(const Vector3f &wo,
                                                Vector2f sample) const {
  auto localToWorld = [](const Vector3f &wo, const Vector3f &local) {
    Vector3f tangent{1.f, 0.f, .0f};
    Vector3f bitangent;
    if (std::abs(dot(tangent, wo)) > .9f) {
      tangent = Vector3f(.0f, 1.f, .0f);
    }
    bitangent = normalize(cross(tangent, wo));
    tangent = normalize(cross(wo, bitangent));

    Vector3f world = local[0] * tangent + local[1] * wo + local[2] * bitangent;
    return world;
  };

  float cosTheta, sinTheta, phi = 2 * PI * sample[1];
  if (std::abs(g) < 1e-3f)
    cosTheta = 1 - 2 * sample[0];
  else {
    cosTheta = -1 / (2 * g) *
               (1 + g * g -
                ((1 - g * g) / (1 + g - 2 * g * sample[0])) *
                    ((1 - g * g) / (1 + g - 2 * g * sample[0])));
  }
  sinTheta = fm::sqrt(std::max(.0f, 1 - cosTheta * cosTheta));
  Vector3f local{sinTheta * fm::sin(phi), cosTheta, sinTheta * fm::cos(phi)};
  Vector3f world = localToWorld(wo, local);
  return {Spectrum(1.f), world, pdf(world, wo), false};
}