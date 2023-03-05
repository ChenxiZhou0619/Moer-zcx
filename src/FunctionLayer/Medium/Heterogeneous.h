#include "Medium.h"
#include <openvdb/openvdb.h>
class HeterogeneousMedium : public Medium {
public:
  HeterogeneousMedium() = delete;

  HeterogeneousMedium(const Json &json);

  virtual Spectrum sample(const Ray &ray, float tmax, Vector2f sample,
                          MediumIntersection *mits) const override;

  virtual Spectrum Tr(Point3f origin, Vector3f direction,
                      float distance) const override;

  virtual Spectrum SigmaS(Point3f position) const override;

public:
  Point3f boxMin, boxMax;

private:
  openvdb::FloatGrid::Ptr density;
};