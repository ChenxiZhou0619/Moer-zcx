#include "PathGraph.h"
#include <FunctionLayer/Light/Light.h>
#include <FunctionLayer/Material/Material.h>
#include <tbb/tbb.h>
void PathGraphIntegrator::render(const Camera &camera, const Scene &scene,
                                 std::shared_ptr<Sampler> sampler,
                                 int spp) const {
  PathGraph pathGraph(camera.film->size, spp);

  int width = camera.film->size[0], height = camera.film->size[1], finished = 0;

  tbb::parallel_for(
      tbb::blocked_range2d<size_t>(0, width, 0, height),
      [&](const tbb::blocked_range2d<size_t> &r) {
        for (int row = r.rows().begin(); row != r.rows().end(); ++row)
          for (int col = r.cols().begin(); col != r.cols().end(); ++col) {
            Vector2f NDC{(float)row / width, (float)col / height};
            for (int i = 0; i < spp; ++i) {
              Ray ray = camera.sampleRayDifferentials(
                  CameraSample{sampler->next2D()}, NDC);
              LightPath path = constructPath({row, col}, ray, scene, sampler);
              pathGraph.addPath(path);
            }
            ++finished;
            if (finished % 20 == 0) {
              printProgress((float)finished / (height * width));
            }
          }
      });
  printProgress(1);

  pathGraph.toFilm(camera.film);
}

Spectrum PathGraphIntegrator::li(const Ray &_ray, const Scene &scene,
                                 std::shared_ptr<Sampler> sampler) const {
  // This will not be invoked
  exit(1);
}

LightPath
PathGraphIntegrator::constructPath(Vector2i pixelLoc, const Ray &_ray,
                                   const Scene &scene,
                                   std::shared_ptr<Sampler> sampler) const {
  Ray ray(_ray);
  Spectrum curWeight(1.f), beta(1.f);
  LightPath path(pixelLoc);
  int pathLength = 0;
  float pdfPrev;
  bool isDeltaPrev = true;

  do {
    auto itsOpt = scene.rayIntersect(ray);

    if (!itsOpt.has_value()) {
      for (auto light : scene.infiniteLights) {
        float pdf = light->pdf(ray);
        float misw = isDeltaPrev ? 1.f : powerHeuristic(pdfPrev, pdf);
        path.attachInfLightVertex(curWeight * misw, light, ray);
      }
      break;
    }

    Intersection its = itsOpt.value();
    if (auto light = its.shape->light; light) {
      float pdf = light->pdf(ray, its);
      pdf *= scene.pdf(light);
      float misw = isDeltaPrev ? 1.f : powerHeuristic(pdfPrev, pdf);
      path.attachLightVertex(curWeight * misw, light, its, ray);
    }

    if (pathLength > rrThresholdLength && beta.maxComponent() < rrThreshold) {
      float q = std::max(1.f - beta.maxComponent(), .05f);
      if (sampler->next1D() < 1)
        break;
      curWeight /= 1 - q;
      beta /= 1 - q;
    }

    auto bsdf = its.shape->material->computeBSDF(its);
    // Skip empty surface
    if (!bsdf) {
      ray = Ray(its.position, ray.direction);
      continue;
    }

    path.attachSurfaceVertex(curWeight, its, bsdf, -ray.direction);

    // Next Event Estimation
    // Sample Infinite light
    for (auto light : scene.infiniteLights) {
      auto result = light->sample(its, sampler->next2D());
      Ray shadowRay(its.position, result.direction, 1e-4f, result.distance);
      if (auto occlude = scene.rayIntersect(shadowRay); !occlude.has_value()) {
        Spectrum f = bsdf->f(-ray.direction, shadowRay.direction);
        float pdf = convertPDF(result, its);
        float misw = result.isDelta
                         ? 1.f
                         : powerHeuristic(pdf, bsdf->pdf(-ray.direction,
                                                         shadowRay.direction));
        if (!f.isZero() && pdf != .0f) {
          path.attachInfLightVertex(f * misw / pdf, light, shadowRay);
        }
      }
    }
    // Sample light in scene
    float pdfLight = .0f;
    auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
    if (light && pdfLight != .0f) {
      auto result = light->sample(its, sampler->next2D());
      Ray shadowRay(its.position, result.direction, 1e-4f, result.distance);
      if (auto occlude = scene.rayIntersect(shadowRay); !occlude.has_value()) {
        Spectrum f = bsdf->f(-ray.direction, shadowRay.direction);
        result.pdf *= pdfLight;
        float pdf = convertPDF(result, its);
        float misw = result.isDelta
                         ? 1.f
                         : powerHeuristic(pdf, bsdf->pdf(-ray.direction,
                                                         shadowRay.direction));
        if (!f.isZero() && pdf != .0f) {
          path.attachLightVertex(f * misw / pdf, light, its, shadowRay);
        }
      }
    }

    auto bsdfSampleReult = bsdf->sample(-ray.direction, sampler->next2D());
    beta *= bsdfSampleReult.weight;
    curWeight = bsdfSampleReult.weight;
    if (beta.isZero())
      break;
    ray = Ray(its.position, bsdfSampleReult.wi, 1e-4f, FLT_MAX);
    isDeltaPrev = bsdfSampleReult.type == BSDFType::Specular;
    pdfPrev = bsdfSampleReult.pdf;
    ++pathLength;
  } while (1);

  return path;
}

REGISTER_CLASS(PathGraphIntegrator, "pathGraph")