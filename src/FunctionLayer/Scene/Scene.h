#pragma once
#include <FunctionLayer/Acceleration/Acceleration.h>
#include <FunctionLayer/Light/EnvironmentLight.h>
#include <FunctionLayer/Light/Light.h>
#include <FunctionLayer/Medium/Medium.h>
#include <ResourceLayer/JsonUtil.h>
class Scene {
public:
  Scene() = delete;

  Scene(const Json &json);

  ~Scene() = default;

  std::optional<Intersection> rayIntersect(const Ray &ray) const;

  std::shared_ptr<Light> sampleLight(float sample, float *pdf) const;

  float pdf(std::shared_ptr<Light> light) const;

  std::vector<std::shared_ptr<InfiniteLight>> infiniteLights;

  std::shared_ptr<Medium> sceneMedium;

private:
  std::shared_ptr<Acceleration> acceleration;

  Distribution1D<std::shared_ptr<Light>> lightDistribution;
};
