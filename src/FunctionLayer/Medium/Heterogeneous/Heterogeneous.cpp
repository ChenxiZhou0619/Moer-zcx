#include "Heterogeneous.h"
#include "RegularTracker.h"
#include <FunctionLayer/Shape/Intersection.h>
#include <ResourceLayer/FileUtil.h>
#include <nanovdb/util/IO.h>
#include <nanovdb/util/SampleFromVoxels.h>

void DEBUG(nanovdb::Vec3<float> v) {
  printf("[%.5f, %.5f, %.5f]\n", v[0], v[1], v[2]);
}

void DEBUG(nanovdb::Vec3<double> v) {
  printf("[%.5f, %.5f, %.5f]\n", v[0], v[1], v[2]);
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

  auto voxelSize = densityFloatGrid->voxelSize();
  if (voxelSize[0] != voxelSize[1] || voxelSize[0] != voxelSize[2]) {
    std::cerr << "Only support cube voxel!\n";
    exit(1);
  }
  voxelScale = voxelSize[0];

  minIndex[0] = densityFloatGrid->indexBBox().min().x();
  minIndex[1] = densityFloatGrid->indexBBox().min().y();
  minIndex[2] = densityFloatGrid->indexBBox().min().z();

  maxIndex[0] = densityFloatGrid->indexBBox().max().x();
  maxIndex[1] = densityFloatGrid->indexBBox().max().y();
  maxIndex[2] = densityFloatGrid->indexBBox().max().z();
}

//* DDA-like tracking

bool HeterogeneousMedium::Sample_RegularTracking(Ray ray, float tmax,
                                                 Vector2f sample,
                                                 MediumIntersection *mits,
                                                 Spectrum *Tr,
                                                 float *pdf) const {
  Point3f o = ray.origin;
  Vector3f d = ray.direction;

  auto o_grid =
           densityFloatGrid->worldToIndexF(nanovdb::Vec3f(o[0], o[1], o[2])),
       d_gird =
           densityFloatGrid->worldToIndexDirF(nanovdb::Vec3f(d[0], d[1], d[2]));

  o = Point3f{o_grid[0], o_grid[1], o_grid[2]};
  d = normalize(Vector3f{d_gird[0], d_gird[1], d_gird[2]});

  // * thick = \sum sigma_t * step
  float thick = -std::log(1 - sample[0]);

  int index[3];
  float dt, sum = .0f;
  float t_world = .0f;

  //* Init regularTracker (all params in indexSpace)
  RegularTracker rt(minIndex, maxIndex, o, d, tmax / voxelScale);

  // index and dt is in indexSpace
  while (rt.track(index, &dt)) {
    nanovdb::Vec3<float> voxel_loc(index[0] + .5f, index[1] + .5f,
                                   index[2] + .5f);

    float density = GridSampler(densityFloatGrid->tree())(voxel_loc),
          delta = density * dt * voxelScale;

    if (sum + delta >= thick) {

      //* Sample a valid point before exit the medium
      dt = (thick - sum) / density / voxelScale;

      t_world += dt * voxelScale;

      mits->position = ray.at(t_world);

      auto dist = (mits->position - ray.origin).length();

      mits->mp.phase = phase;
      mits->mp.sigma_a = Spectrum(.0f); // TODO
      mits->mp.sigma_s = density;

      *Tr = Spectrum(std::exp(-thick));
      *pdf = (*Tr)[0] * density;

      return true;
    }

    sum += delta;

    t_world += dt * voxelScale;
  }

  *Tr = Spectrum(std::exp(-sum));
  *pdf = (*Tr)[0];

  return false;
}

Spectrum HeterogeneousMedium::Transmittance_RegularTracking(Ray ray,
                                                            float t) const {

  Point3f o = ray.origin;
  Vector3f d = ray.direction;

  auto o_grid =
           densityFloatGrid->worldToIndexF(nanovdb::Vec3f(o[0], o[1], o[2])),
       d_gird =
           densityFloatGrid->worldToIndexDirF(nanovdb::Vec3f(d[0], d[1], d[2]));

  o = Point3f{o_grid[0], o_grid[1], o_grid[2]};
  d = normalize(Vector3f{d_gird[0], d_gird[1], d_gird[2]});

  RegularTracker rt(minIndex, maxIndex, o, d, t / voxelScale);
  int index[3];
  float dt, thick = .0f;

  while (rt.track(index, &dt)) {
    nanovdb::Vec3<float> voxel_loc(index[0] + .5f, index[1] + .5f,
                                   index[2] + .5f);
    float density = GridSampler(densityFloatGrid->tree())(voxel_loc);
    thick += density * dt * voxelScale;
  }
  return Spectrum(std::exp(-thick));
}

REGISTER_CLASS(HeterogeneousMedium, "heterogeneous")