#pragma once
#include <CoreLayer/ColorSpace/Spectrum.h>
#include <FunctionLayer/Camera/Camera.h>
#include <FunctionLayer/Ray/Ray.h>
#include <FunctionLayer/Sampler/Sampler.h>
#include <FunctionLayer/Scene/Scene.h>
#include <ResourceLayer/Factory.h>
#include <ResourceLayer/JsonUtil.h>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

inline void printProgress(float percentage) {
  int val = (int)(percentage * 100);
  int lpad = (int)(percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}

class Integrator {
public:
  Integrator() = default;

  virtual ~Integrator() = default;

  Integrator(const Json &json) {}

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const = 0;

  virtual void render(const Camera &camera, const Scene &scene,
                      std::shared_ptr<Sampler> sampler, int spp) const = 0;
  ;
};

class PixelIntegrator : public Integrator {
public:
  PixelIntegrator() = default;

  virtual ~PixelIntegrator() = default;

  PixelIntegrator(const Json &json) : Integrator(json) {}

  virtual Spectrum li(const Ray &ray, const Scene &scene,
                      std::shared_ptr<Sampler> sampler) const override = 0;

  virtual void render(const Camera &camera, const Scene &scene,
                      std::shared_ptr<Sampler> sampler, int spp) const override;
};

//* 将result中各种测度下的pdf都转化为立体角测度下的pdf
inline float convertPDF(const LightSampleResult &result,
                        const Intersection &intersection) {
  float pdf = result.pdf;
  float disance = result.distance;
  switch (result.type) {
  case LightType::SpotLight:
    pdf *= disance * disance;
    break;
  case LightType::AreaLight:
    pdf *= disance * disance;
    pdf /= std::abs(dot(result.normal, result.direction));
    break;
  //* 环境光的pdf转换在采样时已经完成
  case LightType::EnvironmentLight:
    break;
  }
  return pdf;
}

inline float powerHeuristic(float pdfA, float pdfB) {
  pdfA *= pdfA;
  pdfB *= pdfB;
  return pdfA / (pdfA + pdfB);
}