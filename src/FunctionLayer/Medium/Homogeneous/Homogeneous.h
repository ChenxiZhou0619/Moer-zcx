#include "../Medium.h"

class HomogeneousMedium : public Medium {
public:
  HomogeneousMedium() = delete;

  HomogeneousMedium(const Json &json);

  virtual bool Sample_RegularTracking(Ray ray, float tmax, Vector2f sample,
                                      MediumIntersection *mits, Spectrum *Tr,
                                      float *pdf) const override;

  virtual Spectrum Transmittance_RegularTracking(Ray ray,
                                                 float t) const override;

private:
  Spectrum sigmaA, sigmaS;
  float sigmaT; // TODO
};