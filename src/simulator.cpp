#include <FunctionLayer/Acceleration/AABB.h>
#include <FunctionLayer/Acceleration/EmbreeBVH.h>
#include <FunctionLayer/Shape/Sphere.h>
#include <ResourceLayer/Factory.h>
#include <ResourceLayer/FileUtil.h>
#include <ResourceLayer/JsonUtil.h>

#include <FunctionLayer/Sampler/IndependentSampler.h>

#include <fstream>
#include <tbb/tbb.h>

/**
 * @brief
 *
 * @param v v should be an normalized vector
 * @return std::pair<float, float> theta, phi
 */
std::pair<float, float> toSphericalCoord(Vector3f v) {
  // x = sin(theta) * cos(phi)
  // y = sin(theta) * sin(phi)
  // z = cos(theta)
}

int main(int argc, char **argv) {
  const std::string sceneDir = std::string(argv[1]);
  FileUtil::setWorkingDirectory(sceneDir);
  std::string sceneJsonPath = FileUtil::getFullPath("scene.json");
  std::ifstream fstm(sceneJsonPath);
  Json json = Json::parse(fstm);

  auto acceleration = std::make_shared<EmbreeBVH>();
  auto shapes = json["shapes"];
  AABB bounding;

  const Shape *bound = nullptr;

  for (int i = 0; i < shapes.size(); ++i) {
    auto shape = Factory::construct_class<Shape>(shapes[i]);
    acceleration->attachShape(shape);
    auto [min, max] = shape->getAABB();
    bounding.Expand(min);
    bounding.Expand(max);

    if (fetchOptional<bool>(shapes[i], "isBound", false)) {
      bound = shape.get();
    }
  }

  //* 构建加速结构
  acceleration->build();

  auto sampler = std::make_shared<IndependentSampler>();

  printf("theta, phi, u, v, val\n");
}