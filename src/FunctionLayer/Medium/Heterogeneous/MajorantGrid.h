#pragma once
#include <CoreLayer/Math/Math.h>
#include <FunctionLayer/Acceleration/AABB.h>

struct MajorantTracker {
public:
  bool track(int index[3], float *dt);

public:
  bool terminate = false;
  float tmin, tmax;
  float nextCrossingT[3], deltaT[3];
  int step[3], voxelLimit[3], voxel[3];

  int debug = 0;
};

struct MajorantGrid {
public:
  MajorantGrid() = default;

  MajorantGrid(int _resolution[3], Point3f boxMin, Point3f boxMax,
               float voxelScale)
      : box(boxMin, boxMax), voxelScale(voxelScale) {
    resulotion[0] = _resolution[0];
    resulotion[1] = _resolution[1];
    resulotion[2] = _resolution[2];
    majorantVoxel =
        std::vector<float>(resulotion[0] * resulotion[1] * resulotion[2]);
    voxelSize[0] = (boxMax[0] - boxMin[0]) / resulotion[0];
    voxelSize[1] = (boxMax[1] - boxMin[1]) / resulotion[1];
    voxelSize[2] = (boxMax[2] - boxMin[2]) / resulotion[2];
  }

  float at(int x, int y, int z) const {
    int offset = x + y * resulotion[0] + z * resulotion[1] * resulotion[0];
    return majorantVoxel[offset];
  }

  void set(int x, int y, int z, float val) {
    int offset = x + y * resulotion[0] + z * resulotion[1] * resulotion[0];
    majorantVoxel[offset] = val;
  }

  MajorantTracker getTracker(Point3f origin_u, Vector3f dir_u,
                             float t_grid) const;

public:
  AABB box;
  int resulotion[3];
  float voxelSize[3], voxelScale;
  std::vector<float> majorantVoxel;
};