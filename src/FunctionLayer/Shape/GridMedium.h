#pragma once
#include "Cube.h"
#include <FunctionLayer/Medium/Heterogeneous.h>
class GridMedium : public Cube {
public:
  GridMedium() = delete;

  GridMedium(const Json &json);

protected:
  std::shared_ptr<HeterogeneousMedium> grid;
};