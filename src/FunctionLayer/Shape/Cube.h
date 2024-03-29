#pragma once
#include "Shape.h"

class Cube : public Shape {
public:
  Cube() = delete;

  Cube(const Json &json);

  Cube(Point3f _boxMin, Point3f _boxMax, Transform transform);

  virtual bool rayIntersectShape(const Ray &ray, float *distance, int *primID,
                                 float *u, float *v) const override;

  virtual void fillIntersection(float distance, int primID, float u, float v,
                                Intersection *intersection) const override;

  virtual void uniformSampleOnSurface(Vector2f sample,
                                      Intersection *intersection,
                                      float *pdf) const override {
    // TODO finish this
    return;
  }

protected:
  Point3f boxMin, boxMax;
};