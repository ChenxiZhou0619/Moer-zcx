#pragma once

#include "Geometry.h"
#include <functional>
#include <vector>

template <typename T> class Distribution1D {
public:
  Distribution1D() = default;

  Distribution1D(std::size_t size) {
    data.reserve(size);
    cdf.reserve(size + 1);
    cdf.emplace_back(.0f);
  }

  void add(const T &item, float weight) {
    data.emplace_back(item);
    cdf.emplace_back(weight + cdf.back());
  }

  void build() {
    sum = cdf.back();
    float invTotal = 1.f / sum;
    for (auto &pdf : cdf) {
      pdf *= invTotal;
    }
  }

  T sample(float sample, float *pdf) const {
    if (cdf.size() == 1) {
      *pdf = .0f;
      return T();
    }

    auto entry = std::lower_bound(cdf.cbegin(), cdf.cend(), sample);
    size_t index = entry - cdf.cbegin() - 1;
    *pdf = cdf[index + 1] - cdf[index];
    return data[std::min(index, cdf.size() - 2)];
  }

  float pdf(T sampled) const {
    auto entry = std::find(data.cbegin(), data.cend(), sampled);
    if (entry == data.cend())
      return .0f;
    size_t index = entry - data.cbegin();
    return cdf[index + 1] - cdf[index];
  }

  float sum = .0f;

protected:
  std::vector<T> data;
  std::vector<float> cdf;
};
