#include "Heterogeneous.h"
#include "IsotropicPhase.h"
#include <FunctionLayer/Sampler/IndependentSampler.h>
#include <FunctionLayer/Shape/Intersection.h>
#include <ResourceLayer/FileUtil.h>
#include <openvdb/tools/ValueTransformer.h>

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

  openvdb::tools::foreach (density->beginValueOn(),
                           [&](const openvdb::FloatGrid::ValueOnIter &itr) {
                             itr.setValue(*itr * 15.f);
                             sigmaTMax = sigmaTMax > *itr ? sigmaTMax : *itr;
                           });

  std::cout << sigmaTMax << std::endl;

  auto bound = density->evalActiveVoxelBoundingBox();

  auto min = density->indexToWorld(bound.min()),
       max = density->indexToWorld(bound.max());

  boxMin = Point3f(min.x(), min.y(), min.z());
  boxMax = Point3f(max.x(), max.y(), max.z());

  phase = std::make_shared<IsotropicPhase>();
}

Spectrum HeterogeneousMedium::sample(const Ray &ray, float tmax,
                                     Vector2f sample,
                                     MediumIntersection *mits) const {
  const static openvdb::tools::GridSampler<openvdb::FloatGrid,
                                           openvdb::tools::BoxSampler>
      gridSampler(*density);
  static IndependentSampler sampler;

  if (tmax > 1e5f)
    return Spectrum(1.f);

  float invSigmaTTmax = 1.f / sigmaTMax, distance = .0f;
  Spectrum tr(1.f);
  do {
    distance += -fm::log(1 - sampler.next1D()) * invSigmaTTmax;
    if (distance >= tmax) {
      mits->medium = nullptr;
      mits->valid = false;
      mits->position = ray.at(tmax);
      mits->distance = tmax;
      tr *= Tr(ray.origin, ray.direction, tmax);
      for (int i = 0; i < 3; ++i) {
        mits->pdf += tr[i];
      }
      mits->pdf /= 3.f;
      break;
    }
    Point3f position = ray.at(distance);
    float sigmaT = gridSampler.wsSample(
        openvdb::Vec3R(position[0], position[1], position[2]));
    if (sampler.next1D() < sigmaT * invSigmaTTmax) {
      // Yes, sample a medium scatter
      mits->medium = this;
      mits->valid = true;
      mits->position = position;
      mits->distance = distance;
      mits->sigmaS = sigmaT;
      tr = Tr(ray.origin, ray.direction, distance);
      for (int i = 0; i < 3; ++i) {
        mits->pdf += sigmaT * tr[i];
      }
      mits->pdf /= 3.f;
      break;
    }
  } while (1);
  return tr / mits->pdf;
}

Spectrum HeterogeneousMedium::Tr(Point3f origin, Vector3f direction,
                                 float distance) const {
  const static openvdb::tools::GridSampler<openvdb::FloatGrid,
                                           openvdb::tools::BoxSampler>
      gridSampler(*density);
  static IndependentSampler sampler;
  float invSigmaTTmax = 1.f / sigmaTMax, t = .0f;
  Spectrum tr(1.f);

  int step = 0;
  do {
    t += -fm::log(1 - sampler.next1D()) * invSigmaTTmax;
    if (t > distance)
      break;
    Point3f position = origin + direction * t;
    float sigmaT = gridSampler.wsSample(
        openvdb::Vec3R(position[0], position[1], position[2]));
    tr *= Spectrum(1.f - std::max(.0f, sigmaT * invSigmaTTmax));
    if (tr.isZero())
      break;
  } while (++step < 100);
  return tr;
}

Spectrum HeterogeneousMedium::SigmaS(Point3f position) const {
  const static openvdb::tools::GridSampler<openvdb::FloatGrid,
                                           openvdb::tools::BoxSampler>
      gridSampler(*density);
  return gridSampler.wsSample(
      openvdb::Vec3R(position[0], position[1], position[2]));
}

REGISTER_CLASS(HeterogeneousMedium, "heterogeneous")