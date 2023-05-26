#pragma once
#include "Texture.h"
#include <CoreLayer/ColorSpace/Spectrum.h>

class CheckboardTexture : public Texture<Spectrum> {
public:
  CheckboardTexture(const Spectrum &_data1, const Spectrum &_data2, float _size)
      : data1(_data1), data2(_data2), size(_size) {}

  CheckboardTexture(const Json &json);

  virtual Spectrum evaluate(const Intersection &intersection) const override {
    TextureCoord texCoord{intersection.texCoord, Vector2f{.0f, .0f},
                          Vector2f{.0f, .0f}};
    return evaluate(texCoord);
  }

  virtual Spectrum evaluate(const TextureCoord &texCoord) const override {
    int xIdx = texCoord.coord[0] / size, yIdx = texCoord.coord[1] / size;
    return ((xIdx + yIdx) % 2) ? data1 : data2;
  }

private:
  float size;
  Spectrum data1, data2;
};