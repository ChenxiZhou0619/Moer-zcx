#include "Heterogeneous.h"
#include <ResourceLayer/FileUtil.h>

HeterogeneousMedium::HeterogeneousMedium(const Json &json) : Medium(json) {
  std::string path = fetchRequired<std::string>(json, "file");
  path = FileUtil::getFullPath(path);

  openvdb::initialize();
  openvdb::io::File vdbFile(path);
  vdbFile.open();

  auto itr = vdbFile.beginName();
  if (itr == vdbFile.endName()) {
    std::cout << "Error! No grid in file\n";
    exit(1);
  }

  density = openvdb::gridPtrCast<openvdb::FloatGrid>(
      vdbFile.readGrid(itr.gridName()));

  auto bound = density->evalActiveVoxelBoundingBox();

  auto min = density->indexToWorld(bound.min()),
       max = density->indexToWorld(bound.max());

  boxMin = Point3f(min.x(), min.y(), min.z());
  boxMax = Point3f(max.x(), max.y(), max.z());
}

Spectrum HeterogeneousMedium::sample(const Ray &ray, float tmax,
                                     Vector2f sample,
                                     MediumIntersection *mits) const {
  //
}

Spectrum HeterogeneousMedium::Tr(Point3f origin, Vector3f direction,
                                 float distance) const {
  //
}

Spectrum HeterogeneousMedium::SigmaS(Point3f position) const {
  //
}

REGISTER_CLASS(HeterogeneousMedium, "heterogeneous")