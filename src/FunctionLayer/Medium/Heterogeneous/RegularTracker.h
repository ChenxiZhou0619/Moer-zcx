#pragma once
#include <CoreLayer/Math/Math.h>
#include <FunctionLayer/Ray/Ray.h>
class RegularTracker {
public:
  RegularTracker() = delete;

  RegularTracker(const int min[3], const int max[3], Point3f origin,
                 Vector3f direction, float tmax, float voxelScale);

  bool track(int index[3], float *dt);

public:
  bool terminate = false;
  float voxelScale;
  float tmin, tmax;
  float nextCrossingT[3], deltaT[3];
  int step[3], voxelLimit[3], voxel[3];
};