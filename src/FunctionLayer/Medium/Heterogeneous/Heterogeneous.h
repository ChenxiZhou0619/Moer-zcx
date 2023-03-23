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

public:
  Point3f boxMin, boxMax; // Bound the volume

private:
  nanovdb::GridHandle<BufferT> densityGrid;
  nanovdb::GridHandle<BufferT> temperatureGrid;

  const nanovdb::FloatGrid *densityFloatGrid = nullptr;
  const nanovdb::FloatGrid *temperatureFloatGrid = nullptr;

  float voxelScale;
  int minIndex[3], maxIndex[3];

  float sigma_maj = .0f; // TODO replace with maj_grid
  MajorantGrid majorantGrid;
};