#pragma once
#include <CoreLayer/Math/Math.h>
//* 使用RGB三通道颜色空间

class SpectrumRGB {
public:
  SpectrumRGB() : rgb(0.f) {}

  SpectrumRGB(float f) : rgb(f) {}

  SpectrumRGB(float r, float g, float b) : rgb(r, g, b) {}

  SpectrumRGB(Vector3f _rgb) : rgb(_rgb) {}

  SpectrumRGB operator+(const SpectrumRGB &rhs) const {
    return SpectrumRGB(rgb + rhs.rgb);
  }

  SpectrumRGB &operator+=(const SpectrumRGB &rhs) {
    rgb += rhs.rgb;
    return *this;
  }

  SpectrumRGB operator-(const SpectrumRGB &rhs) const {
    return SpectrumRGB(rgb - rhs.rgb);
  }

  SpectrumRGB &operator-=(const SpectrumRGB &rhs) {
    rgb -= rhs.rgb;
    return *this;
  }

  SpectrumRGB operator*(const SpectrumRGB &rhs) const {
    return SpectrumRGB(rgb * rhs.rgb);
  }

  SpectrumRGB &operator*=(const SpectrumRGB &rhs) {
    rgb *= rhs.rgb;
    return *this;
  }

  SpectrumRGB operator*(float f) const { return SpectrumRGB(rgb * f); }

  SpectrumRGB &operator*=(float f) {
    rgb *= f;
    return *this;
  }

  SpectrumRGB operator/(const SpectrumRGB &rhs) const {
    return SpectrumRGB(rgb / rhs.rgb);
  }

  SpectrumRGB &operator/=(const SpectrumRGB &rhs) {
    rgb /= rhs.rgb;
    return *this;
  }

  SpectrumRGB operator/(float f) const { return SpectrumRGB(rgb / f); }

  SpectrumRGB &operator/=(float f) {
    rgb /= f;
    return *this;
  }

  float operator[](int i) const { return rgb[i]; }

  float &operator[](int i) { return rgb[i]; }

  bool isZero() const { return rgb.isZero(); }

  bool hasNaN() const {
    return std::isnan(rgb[0]) || std::isnan(rgb[1]) || std::isnan(rgb[2]);
  }

  bool hasInf() const {
    return std::isinf(rgb[0]) || std::isinf(rgb[1]) || std::isinf(rgb[2]);
  }

  float maxComponent() const {
    return std::max(std::max(rgb[0], rgb[1]), rgb[2]);
  }

  void debugPrint() const {
    printf("[rgb](");
    for (int i = 0; i < 3; ++i) {
      std::cout << (i == 0 ? '\0' : ',') << rgb[i];
    }
    printf(")%c", '\n');
    fflush(stdout);
  }

private:
  Vector3f rgb;
};

inline SpectrumRGB operator*(float f, const SpectrumRGB &spectrum) {
  return spectrum * f;
}

inline Vector3f toVec3(const SpectrumRGB &spectrum) {
  return Vector3f{spectrum[0], spectrum[1], spectrum[2]};
}

inline SpectrumRGB toSpectrum(const Vector3f &vec) {
  return SpectrumRGB(vec[0], vec[1], vec[2]);
}