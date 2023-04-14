#include "Heterogeneous.h"
#include "RegularTracker.h"
#include <FunctionLayer/Sampler/IndependentSampler.h>
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

float HeterogeneousMedium::scaleSample(nanovdb::Vec3R index,
                                       const nanovdb::FloatGrid *grid) const {

  return GridSampler(grid->tree())(index) * sigmaScale;
}

HeterogeneousMedium::HeterogeneousMedium(const Json &json) : Medium(json) {
  sigmaScale = fetchOptional(json, "sigmaScale", 1.f);
  albedo = fetchOptional(json, "albedo", 1.f);

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

  //* Init majorant grid
  int resolution[] = {32, 32, 32};
  majorantGrid = MajorantGrid(resolution, boxMin, boxMax, voxelScale);

  constexpr float delta = 1.f;
  //* Set majorant at each majorantVoxel
  for (int iz = 0; iz < resolution[2]; ++iz) {
    for (int iy = 0; iy < resolution[1]; ++iy) {
      for (int ix = 0; ix < resolution[0]; ++ix) {
        Point3f u_min, u_max;
        u_min[0] = (float)ix / resolution[0],
        u_min[1] = (float)iy / resolution[1];
        u_min[2] = (float)iz / resolution[2];

        u_max[0] = (float)(ix + 1.f) / resolution[0],
        u_max[1] = (float)(iy + 1.f) / resolution[1];
        u_max[2] = (float)(iz + 1.f) / resolution[2];

        auto [vmin, vmax] = majorantGrid.box.SubBox(u_min, u_max);
        nanovdb::Vec3<float> voxelMin = densityFloatGrid->worldToIndexF(
            nanovdb::Vec3<float>(vmin[0], vmin[1], vmin[2]));
        nanovdb::Vec3<float> voxelMax = densityFloatGrid->worldToIndexF(
            nanovdb::Vec3<float>(vmax[0], vmax[1], vmax[2]));
        float maj = .0f;
        for (int nz = voxelMin[2] - delta; nz <= voxelMax[2] + delta; ++nz)
          for (int ny = voxelMin[1] - delta; ny <= voxelMax[1] + delta; ++ny)
            for (int nx = voxelMin[0] - delta; nx <= voxelMax[0] + delta; ++nx)
              maj = std::max(maj, scaleSample(nanovdb::Vec3R(nx, ny, nz),
                                              densityFloatGrid));
        majorantGrid.set(ix, iy, iz, maj);
      }
    }
  }
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
  RegularTracker rt(minIndex, maxIndex, o, d, tmax, voxelScale);

  // dt in worldSpace
  while (rt.track(index, &dt)) {
    nanovdb::Vec3<double> voxel_loc(index[0] + .5f, index[1] + .5f,
                                    index[2] + .5f);

    float density = scaleSample(voxel_loc, densityFloatGrid),
          delta = density * dt;

    if (sum + delta >= thick) {

      //* Sample a valid point before exit the medium
      dt = (thick - sum) / density;

      t_world += dt;

      mits->position = ray.at(t_world);
      mits->mp.phase = phase;
      mits->mp.sigma_s = density * albedo;
      mits->mp.sigma_a = Spectrum(density) - mits->mp.sigma_s;

      *Tr = Spectrum(std::exp(-thick));
      *pdf = (*Tr)[0] * density;

      return true;
    }

    sum += delta;
    t_world += dt;
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

  RegularTracker rt(minIndex, maxIndex, o, d, t, voxelScale);
  int index[3];
  float dt, thick = .0f;

  while (rt.track(index, &dt)) {
    nanovdb::Vec3<double> voxel_loc(index[0] + .5f, index[1] + .5f,
                                    index[2] + .5f);
    float density = scaleSample(voxel_loc, densityFloatGrid);
    thick += density * dt;
  }
  return Spectrum(std::exp(-thick));
}

bool HeterogeneousMedium::Sample_MajorantTracking(Ray ray, float tmax,
                                                  Vector2f sample,
                                                  MediumIntersection *mits,
                                                  Spectrum *Tr,
                                                  float *pdf) const {
  *Tr = Spectrum(1.f);
  *pdf = 1.f;

  Point3f u_coord = majorantGrid.box.UniformCoord(ray.origin);

  MajorantTracker mt = majorantGrid.getTracker(u_coord, ray.direction, tmax);

  float thick = -std::log(1 - sample[0]);

  int index[3];
  float dt, sum = .0f;
  float t_world = .0f;

  while (mt.track(index, &dt)) {
    float maj = majorantGrid.at(index[0], index[1], index[2]), delta = maj * dt;

    if (sum + delta >= thick) {
      dt = (thick - sum) / maj;
      t_world += dt;

      mits->position = ray.at(t_world);
      mits->mp.phase = phase;

      nanovdb::Vec3<double> indexLoc =
          densityFloatGrid->worldToIndexF(nanovdb::Vec3<double>(
              mits->position[0], mits->position[1], mits->position[2]));
      float density = scaleSample(indexLoc, densityFloatGrid);
      mits->mp.sigma_maj = maj;
      mits->mp.sigma_s = density * albedo;
      mits->mp.sigma_a = Spectrum(density) - mits->mp.sigma_s;
      return true;
    }

    t_world += dt;
    sum += delta;
  }

  return false;
}

Spectrum HeterogeneousMedium::Transmittance_RatioTracking(Ray ray,
                                                          float t) const {
  static IndependentSampler sampler;

  Point3f u_coord = majorantGrid.box.UniformCoord(ray.origin);
  MajorantTracker mt = majorantGrid.getTracker(u_coord, ray.direction, t);
  Spectrum Tr(1.f);

  int index[3], axis;
  float dt, thick_sum = .0f, t_sum = .0f;

  float thick = -std::log(1 - sampler.next1D());
  while (mt.track(index, &dt, &axis)) {
    float maj = majorantGrid.at(index[0], index[1], index[2]), delta = dt * maj;

    if (thick_sum + delta >= thick) {
      float step = (thick - thick_sum) / maj;
      t_sum += step;
      Point3f position = ray.at(t_sum);
      nanovdb::Vec3<double> indexLoc = densityFloatGrid->worldToIndexF(
          nanovdb::Vec3<double>(position[0], position[1], position[2]));
      float sigma_t = scaleSample(indexLoc, densityFloatGrid);
      Tr *= Spectrum(1.f) - sigma_t / maj;

      // Update
      thick_sum = .0f;
      thick = -std::log(1 - sampler.next1D());
      mt.terminate = false;
      mt.voxel[0] = index[0];
      mt.voxel[1] = index[1];
      mt.voxel[2] = index[2];
      mt.tmin = mt.tmin - dt + step;
      mt.nextCrossingT[axis] -= mt.deltaT[axis];
      continue;
    }

    t_sum += dt;
    thick_sum += delta;
  }

  return Tr;
}

REGISTER_CLASS(HeterogeneousMedium, "heterogeneous")