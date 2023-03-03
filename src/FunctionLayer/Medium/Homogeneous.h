#include "Medium.h"

class HomogeneousMedium : public Medium {
public:
  HomogeneousMedium() = delete;

  HomogeneousMedium(const Json &json);

  virtual Spectrum sample(const Ray &ray, float tmax, Vector2f sample,
                          MediumIntersection *mits) const override;

  virtual Spectrum Tr(Point3f origin, Vector3f direction,
                      float distance) const override;

  virtual Spectrum SigmaS(Point3f position) const override;

private:
  Spectrum sigmaA, sigmaS;
  float sigmaT; // TODO
};