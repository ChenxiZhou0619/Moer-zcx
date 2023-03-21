#include "Heterogeneous.h"
#include "RegularTracker.h"
#include <FunctionLayer/Shape/Intersection.h>
#include <ResourceLayer/FileUtil.h>
#include <nanovdb/util/IO.h>
#include <nanovdb/util/SampleFromVoxels.h>

void DEBUG(nanovdb::Vec3<float> v) {
  printf("[%.2f, %.2f, %.2f]\n", v[0], v[1], v[2]);
}

using GridSampler =
    nanovdb::SampleFromVoxels<nanovdb::FloatGrid::TreeType, 1, false>;

HeterogeneousMedium::HeterogeneousMedium(const Json &json) : Medium(json) {
  std::string filename = fetchRequired<std::string>(json, "file");
  filename = FileUtil::getFullPath(filename);

  densityGrid = nanovdb::io::readGrid(filename, "density", 1);
  densityFloatGrid = densityGrid.grid<float>();
  if (!densityGrid) {
    std::cout << ".nvdb file must contains density grid\n";
    exit(1);
  }

  temperatureGrid = nanovdb::io::readGrid(filename, "temperature", 1);
  temperatureFloatGrid = temperatureGrid.grid<float>();

  //* Compute density grid bound
  auto bbox = densityFloatGrid->worldBBox();
  boxMin = Point3f(bbox.min()[0], bbox.min()[1], bbox.min()[2]);
  boxMax = Point3f(bbox.max()[0], bbox.max()[1], bbox.max()[2]);
}

//* DDA-like tracking

bool HeterogeneousMedium::Sample_RegularTracking(Ray ray, float tmax,
                                                 Vector2f sample,
                                                 MediumIntersection *mits,
                                                 Spectrum *Tr,
                                                 float *pdf) const {
  const static int minIndex[3] = {densityFloatGrid->indexBBox().min().x(),
                                  densityFloatGrid->indexBBox().min().y(),
                                  densityFloatGrid->indexBBox().min().z()},
                   maxIndex[3] = {densityFloatGrid->indexBBox().max().x(),
                                  densityFloatGrid->indexBBox().max().y(),
                                  densityFloatGrid->indexBBox().max().z()};

  Point3f o = ray.origin;
  Vector3f d = ray.direction;

  auto o_grid =
           densityFloatGrid->worldToIndexF(nanovdb::Vec3f(o[0], o[1], o[2])),
       d_gird =
           densityFloatGrid->worldToIndexDirF(nanovdb::Vec3f(d[0], d[1], d[2]));

  float scale = d_gird.length(), invScale = 1.f / scale;

  o = Point3f{o_grid[0], o_grid[1], o_grid[2]};
  d = normalize(Vector3f{d_gird[0], d_gird[1], d_gird[2]});

  // * thick = \sum sigma_t * step
  float thick = -std::log(1 - sample[0]);

  int index[3];
  float dt, sum = .0f;

  //* Init regularTracker (all params in indexSpace)
  RegularTracker rt(minIndex, maxIndex, o, d, tmax * scale);

  // index and dt is in indexSpace
  while (rt.track(index, &dt)) {
    nanovdb::Vec3<float> voxel_loc(index[0] + .5f, index[1] + .5f,
                                   index[2] + .5f);

    float density = GridSampler(densityFloatGrid->tree())(voxel_loc),
          delta = density * dt * invScale;

    if (sum + delta >= thick) {
      //* Sample a valid point before exit the medium
      dt = (thick - sum) / density;

      float p_x = index[0] + d[0] * dt, p_y = index[1] + d[1] * dt,
            p_z = index[2] + d[2] * dt;

      auto p = densityFloatGrid->indexToWorldF(nanovdb::Vec3f(p_x, p_y, p_z));
      mits->position = Point3f{p[0], p[1], p[2]};
      mits->mp.phase = phase;
      mits->mp.sigma_a = Spectrum(.0f); // TODO
      mits->mp.sigma_s = density;

      *Tr = Spectrum(std::exp(-thick));
      *pdf = (*Tr)[0] * density;
      return true;
    }

    sum += delta;
  }

  *Tr = Spectrum(std::exp(-sum));
  *pdf = (*Tr)[0];
  return false;
}

Spectrum HeterogeneousMedium::Transmittance_RegularTracking(Ray ray,
                                                            float t) const {
  const static int minIndex[3] = {densityFloatGrid->indexBBox().min().x(),
                                  densityFloatGrid->indexBBox().min().y(),
                                  densityFloatGrid->indexBBox().min().z()},
                   maxIndex[3] = {densityFloatGrid->indexBBox().max().x(),
                                  densityFloatGrid->indexBBox().max().y(),
                                  densityFloatGrid->indexBBox().max().z()};

  Point3f o = ray.origin;
  Vector3f d = ray.direction;

  auto o_grid =
           densityFloatGrid->worldToIndexF(nanovdb::Vec3f(o[0], o[1], o[2])),
       d_gird =
           densityFloatGrid->worldToIndexDirF(nanovdb::Vec3f(d[0], d[1], d[2]));

  float scale = d_gird.length(), invScale = 1.f / scale;

  o = Point3f{o_grid[0], o_grid[1], o_grid[2]};
  d = normalize(Vector3f{d_gird[0], d_gird[1], d_gird[2]});

  RegularTracker rt(minIndex, maxIndex, o, d, t * scale);
  int index[3];
  float dt, thick = .0f;

  while (rt.track(index, &dt)) {
    nanovdb::Vec3<float> voxel_loc(index[0] + .5f, index[1] + .5f,
                                   index[2] + .5f);
    float density = GridSampler(densityFloatGrid->tree())(voxel_loc);
    thick += density * dt * invScale;
  }
  return Spectrum(std::exp(-thick));
}

REGISTER_CLASS(HeterogeneousMedium, "heterogeneous")