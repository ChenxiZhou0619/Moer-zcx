#include "GridMedium.h"
#include <FunctionLayer/Material/Material.h>

GridMedium::GridMedium(const Json &json) : Cube(json) {
  if (!json.contains("medium") ||
      strcmp(fetchRequired<std::string>(json["medium"], "type").c_str(),
             "heterogeneous") != 0) {
    std::cerr << "Error, gridMedium must contains a heterogeneous medium!\n";
    exit(1);
  }
  grid = Factory::construct_class<HeterogeneousMedium>(json["medium"]);
  this->material->setMedium(grid);

  // reset cube
  boxMin = grid->boxMin;
  boxMax = grid->boxMax;

  // 构建AABB
  pMin = pMax = transform.toWorld(boxMin);
  for (int i = 0; i < 8; ++i) {
    Point3f p;
    p[0] = (i & 0b100) ? boxMax[0] : boxMin[0];
    p[1] = (i & 0b010) ? boxMax[1] : boxMin[1];
    p[2] = (i & 0b001) ? boxMax[2] : boxMin[2];
    p = transform.toWorld(p);

    for (int j = 0; j < 3; ++j) {
      pMin[j] = std::min(pMin[j], p[j]);
      pMax[j] = std::max(pMax[j], p[j]);
    }
  }

  // 在计算时，所有计算都是在局部坐标系内完成的，因此这里只对boxMin和boxMax做scale操作
  Matrix4f scale = transform.scale;
  vecmat::vec4f min{boxMin[0], boxMin[1], boxMin[2], 1.f},
      max{boxMax[0], boxMax[1], boxMax[2], 1.f};
  min = scale * min, max = scale * max;
  min /= min[3], max /= max[3];
  boxMin = Point3f{min[0], min[1], min[2]};
  boxMax = Point3f{max[0], max[1], max[2]};
}

REGISTER_CLASS(GridMedium, "gridMedium")