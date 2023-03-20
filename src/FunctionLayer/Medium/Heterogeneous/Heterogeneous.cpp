#include "Heterogeneous.h"
#include <ResourceLayer/FileUtil.h>
#include <nanovdb/util/IO.h>
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

  auto min = densityFloatGrid->indexBBox().min(),
       max = densityFloatGrid->indexBBox().max();
  printf("index min = [%d, %d, %d]\n", min.x(), min.y(), min.z());
  printf("index max = [%d, %d, %d]\n", max.x(), max.y(), max.z());

  exit(1);
}

bool HeterogeneousMedium::Sample_RegularTracking(Ray ray, float tmax,
                                                 Vector2f sample,
                                                 MediumIntersection *mits,
                                                 Spectrum *Tr,
                                                 float *pdf) const {
  //
}

Spectrum HeterogeneousMedium::Transmittance_RegularTracking(Ray ray,
                                                            float t) const {
  //
}

REGISTER_CLASS(HeterogeneousMedium, "heterogeneous")