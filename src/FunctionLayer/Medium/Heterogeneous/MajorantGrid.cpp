#include "MajorantGrid.h"

bool MajorantTracker::track(int *index, float *dt, int *axis) {

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

  if (axis)
    *axis = stepAxis;

  if (nextCrossingT[stepAxis] > tmax) {
    *dt = tmax - tmin;
    tmin = tmax;
    terminate = true;
  } else {
    *dt = nextCrossingT[stepAxis] - tmin;
    tmin = nextCrossingT[stepAxis];
  }

  index[0] = voxel[0];
  index[1] = voxel[1];
  index[2] = voxel[2];

  voxel[stepAxis] += step[stepAxis];

  if (voxel[stepAxis] == voxelLimit[stepAxis]) {
    terminate = true;
  }

  nextCrossingT[stepAxis] += deltaT[stepAxis];
  return true;
}

MajorantTracker MajorantGrid::getTracker(Point3f origin_u, Vector3f dir_u,
                                         float t_world) const {
  int index[3], step[3], voxelLimit[3];
  float deltaT[3];        // World Space
  float nextCrossingT[3]; // World Space
  MajorantTracker mt;

  mt.tmin = .0f;
  mt.tmax = t_world;

  for (int axis = 0; axis < 3; ++axis) {
    index[axis] =
        clamp<int>(origin_u[axis] * resulotion[axis], 0, resulotion[axis] - 1);
    deltaT[axis] = ((box.pMax[axis] - box.pMin[axis]) / resulotion[axis]) /
                   std::abs(dir_u[axis]);
    if (dir_u[axis] == -.0f)
      dir_u[axis] = .0f;

    if (dir_u[axis] >= 0) {
      nextCrossingT[axis] =
          ((((float)index[axis] + 1.f) / resulotion[axis]) - origin_u[axis]) *
          (box.pMax[axis] - box.pMin[axis]) / dir_u[axis];
      step[axis] = 1;
      voxelLimit[axis] = resulotion[axis];
    } else {
      nextCrossingT[axis] =
          (((float)index[axis] / resulotion[axis]) - origin_u[axis]) *
          (box.pMax[axis] - box.pMin[axis]) / dir_u[axis];
      step[axis] = -1;
      voxelLimit[axis] = -1;
    }

    mt.voxel[axis] = index[axis];
    mt.deltaT[axis] = deltaT[axis];
    mt.nextCrossingT[axis] = nextCrossingT[axis];
    mt.step[axis] = step[axis];
    mt.voxelLimit[axis] = voxelLimit[axis];
  }

  return mt;
}