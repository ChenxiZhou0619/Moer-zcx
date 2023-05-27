#pragma once
#include <CoreLayer/ColorSpace/Spectrum.h>
#include <FunctionLayer/Shape/Intersection.h>

enum class BSDFType { Diffuse, Specular };

struct BSDFSampleResult {
  Spectrum weight;
  Vector3f wi;
  float pdf;
  BSDFType type;
};

class BSDF {
public:
  BSDF(const Vector3f &_normal, const Vector3f &_tangent,
       const Vector3f &_bitangent, BSDFType _type) {
    normal = _normal;
    tangent = _tangent;
    bitangent = _bitangent;
    type = _type;
  }
  virtual Spectrum f(const Vector3f &wo, const Vector3f &wi) const = 0;
  virtual BSDFSampleResult sample(const Vector3f &wo,
                                  const Vector2f &sample) const = 0;
  virtual float pdf(const Vector3f &wo, const Vector3f &wi) const = 0;

public:
  Vector3f normal, tangent, bitangent; // 构成局部坐标系

  BSDFType type; // bsdf类型

protected:
  Vector3f toLocal(Vector3f world) const {
    return Vector3f{dot(tangent, world), dot(normal, world),
                    dot(bitangent, world)};
  }
  Vector3f toWorld(Vector3f local) const {
    return local[0] * tangent + local[1] * normal + local[2] * bitangent;
  }
};

inline Vector3f Reflect(Vector3f w, Vector3f normal) {
  // assume v is the reflected normal
  // v + w = 2 * cos_v * normal
  float cos = dot(w, normal);
  return normalize(2 * cos * normal - w);
}