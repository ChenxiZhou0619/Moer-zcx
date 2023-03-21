#include "RegularTracker.h"

RegularTracker::RegularTracker(const int min[3], const int max[3],
                               Point3f origin, Vector3f direction,
                               float _tmax) {
  tmin = .0f;
  tmax = _tmax;

  for (int axis = 0; axis < 3; ++axis) {
    voxel[axis] = clamp<int>(origin[axis], min[axis], max[axis]);
    deltaT[axis] = 1.f / std::abs(direction[axis]);

    if (direction[axis] == -.0f)
      direction[axis] = .0f;

    if (direction[axis] >= 0) {
      nextCrossingT[axis] =
          (voxel[axis] + 1.f - origin[axis]) / direction[axis];
      step[axis] = 1;
      voxelLimit[axis] = max[axis] + 1;
    } else {
      nextCrossingT[axis] = (voxel[axis] - origin[axis]) / direction[axis];
      step[axis] = -1;
      voxelLimit[axis] = min[axis] - 1;
    }
  }
}

bool RegularTracker::track(int *index, float *dt) {
  if (terminate)
    return false;

  int stepAxis = -1;

  //* If x hit first
  if (nextCrossingT[0] < nextCrossingT[1] &&
      nextCrossingT[0] < nextCrossingT[2])
    stepAxis = 0;
  //* If y hit first
  else if (nextCrossingT[1] < nextCrossingT[2])
    stepAxis = 1;
  else
    stepAxis = 2;

  if (nextCrossingT[stepAxis] > tmax)
    return false;

  index[0] = voxel[0];
  index[1] = voxel[1];
  index[2] = voxel[2];
  *dt = deltaT[stepAxis];

  voxel[stepAxis] += step[stepAxis];

  if (voxel[stepAxis] == voxelLimit[stepAxis]) {
    terminate = true;
  }

  nextCrossingT[stepAxis] += deltaT[stepAxis];
  return true;
}