#pragma once
#include <CoreLayer/ColorSpace/Spectrum.h>
#include <FunctionLayer/Shape/Intersection.h>

class BSSRDF {
public:
  BSSRDF() {}

  virtual ~BSSRDF() {}

  virtual Spectrum sample_Pi(const Shape *shape, Vector2f u2, float u1,
                             Intersection *pi, float *pdf) const = 0;

protected:
};