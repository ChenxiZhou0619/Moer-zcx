#pragma once
#include "Light.h"
#include <CoreLayer/Math/Distribution.h>
#include <FunctionLayer/Ray/Ray.h>
#include <FunctionLayer/Texture/ImageTexture.h>
class Distribution2D {
public:
  Distribution2D() = default;

  Distribution2D(int sizex, int sizey);

  void addPdf(int dimx, int dimy, float pdf);

  void build();

  Vector2i sample(Vector2f sample, float *pdf) const;

  float pdf(Vector2i index) const;

protected:
  std::vector<float> marginCdf;
  std::vector<std::vector<float>> cdf;
};

inline Distribution2D::Distribution2D(int sizex, int sizey) {
  marginCdf.reserve(sizex + 1);
  marginCdf.emplace_back(.0f);

  for (int i = 0; i < sizex; ++i) {
    std::vector<float> _cdf(sizey + 1);
    _cdf[0] = .0f;
    cdf.emplace_back(_cdf);
  }
}

inline void Distribution2D::addPdf(int dimx, int dimy, float pdf) {
  cdf[dimx][dimy + 1] = cdf[dimx][dimy] + pdf;
}

inline void Distribution2D::build() {
  for (int i = 0; i < cdf.size(); ++i) {
    marginCdf.emplace_back(marginCdf.back() + cdf[i].back());
    float invSum = 1.f / cdf[i].back();
    for (int j = 0; j < cdf[i].size(); ++j) {
      cdf[i][j] *= invSum;
    }
    cdf[i].back() = 1.f;
  }
  float invSum = 1.f / marginCdf.back();
  for (int i = 0; i < marginCdf.size(); ++i) {
    marginCdf[i] *= invSum;
  }
  marginCdf.back() = 1.f;
}

inline Vector2i Distribution2D::sample(Vector2f sample, float *pdf) const {
  float pdf1 = .0f, pdf2 = .0f;
  auto entryx =
      std::lower_bound(marginCdf.cbegin(), marginCdf.cend(), sample[0]);
  int dimx = entryx - marginCdf.cbegin() - 1;
  dimx = clamp(dimx, 0, (int)cdf.size() - 1);
  pdf1 = marginCdf[dimx + 1] - marginCdf[dimx];

  auto entryy =
      std::lower_bound(cdf[dimx].cbegin(), cdf[dimx].cend(), sample[1]);
  int dimy = entryy - cdf[dimx].cbegin() - 1;
  dimy = clamp(dimy, 0, (int)cdf[dimx].size() - 1);
  pdf2 = cdf[dimx][dimy + 1] - cdf[dimx][dimy];
  *pdf = pdf1 * pdf2;
  return {dimx, dimy};
}

inline float Distribution2D::pdf(Vector2i index) const {
  int x = std::min(index[0], (int)cdf.size() - 1),
      y = std::min(index[1], (int)cdf[x].size() - 1);
  float pdfx = marginCdf[x + 1] - marginCdf[x];
  float pdfy = cdf[x][y + 1] - cdf[x][y];
  return pdfx * pdfy;
}

class EnvironmentLight : public InfiniteLight {
public:
  EnvironmentLight() = delete;

  EnvironmentLight(const Json &json);

  virtual Spectrum evaluateEmission(const Ray &ray) const override;

  virtual float pdf(const Ray &ray) const override;

  virtual float pdf(const Ray &ray,
                    const Intersection &its) const override final {
    // This will not be invoke
    return .0f;
  }

  virtual LightSampleResult sample(const Intersection &shadingPoint,
                                   const Vector2f &sample) const override;

private:
  std::shared_ptr<Texture<Spectrum>> environmentMap;
  //  Distribution1D<Vector2i> energyDistribution;
  std::shared_ptr<Distribution2D> energyDist;
};