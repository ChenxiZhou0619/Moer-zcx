#pragma once

#include "../Medium.h"
#include "MajorantGrid.h"
#include <nanovdb/NanoVDB.h>
#include <nanovdb/util/GridHandle.h>

using BufferT = nanovdb::HostBuffer;

class HeterogeneousMedium : public Medium {
public:
  HeterogeneousMedium() = delete;

  HeterogeneousMedium(const Json &json);

  virtual bool Sample_RegularTracking(Ray ray, float tmax, Vector2f sample,
                                      MediumIntersection *mits, Spectrum *Tr,
                                      float *pdf) const override;

  virtual Spectrum Transmittance_RegularTracking(Ray ray,
                                                 float t) const override;

  virtual bool Sample_MajorantTracking(Ray ray, float tmax, Vector2f sample,
                                       MediumIntersection *mits, Spectrum *Tr,
                                       float *pdf) const override;

  virtual Spectrum Transmittance_RatioTracking(Ray ray, float t) const override;

public:
  Point3f boxMin, boxMax; // Bound the volume

private:
  float scaleSample(nanovdb::Vec3R index, const nanovdb::FloatGrid *grid) const;

  nanovdb::GridHandle<BufferT> densityGrid;
  nanovdb::GridHandle<BufferT> temperatureGrid;

  const nanovdb::FloatGrid *densityFloatGrid = nullptr;
  const nanovdb::FloatGrid *temperatureFloatGrid = nullptr;

  float voxelScale;
  int minIndex[3], maxIndex[3];
  MajorantGrid majorantGrid;

  float sigmaScale;
  float albedo;
};