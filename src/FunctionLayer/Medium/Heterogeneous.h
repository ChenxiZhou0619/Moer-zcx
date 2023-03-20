#include "Medium.h"
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
class HeterogeneousMedium : public Medium {
public:
  HeterogeneousMedium() = delete;

  HeterogeneousMedium(const Json &json);

  virtual bool Sample_RegularTracking(Ray ray, float tmax, Vector2f sample,
                                      MediumIntersection *mits, Spectrum *Tr,
                                      float *pdf) const override{
      //
  };

  virtual Spectrum Transmittance_RegularTracking(Ray ray,
                                                 float t) const override{
      //
  };

public:
  Point3f boxMin, boxMax;

private:
  openvdb::FloatGrid::Ptr density;
  float sigmaTMax = .0f;
};