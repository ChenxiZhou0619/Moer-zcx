#include "Integrator.h"

#include <tbb/tbb.h>

#define MULTITHREAD

void PixelIntegrator::render(const Camera &camera, const Scene &scene,
                             std::shared_ptr<Sampler> sampler, int spp) const {
  int width = camera.film->size[0], height = camera.film->size[1], finished = 0;
#ifdef MULTITHREAD
  tbb::parallel_for(
      tbb::blocked_range2d<size_t>(0, width, 0, height),
      [&](const tbb::blocked_range2d<size_t> &r) {
        for (int row = r.rows().begin(); row != r.rows().end(); ++row)
          for (int col = r.cols().begin(); col != r.cols().end(); ++col) {
            Vector2f NDC{(float)row / width, (float)col / height};
            for (int i = 0; i < spp; ++i) {
              Ray ray = camera.sampleRayDifferentials(
                  CameraSample{sampler->next2D()}, NDC);
              Spectrum res = li(ray, scene, sampler);
              if (!res.hasNaN() && !res.hasInf())
                camera.film->deposit({row, col}, res, 1.f);
            }

            ++finished;
            if (finished % 20 == 0) {
              printProgress((float)finished / (height * width));
            }
          }
      });
#else
  for (int row = 0; row != width; ++row)
    for (int col = 0; col != height; ++col) {
      Vector2f NDC{(float)row / width, (float)col / height};
      for (int i = 0; i < spp; ++i) {
        Ray ray =
            camera.sampleRayDifferentials(CameraSample{sampler->next2D()}, NDC);
        camera.film->deposit({row, col}, li(ray, scene, sampler), 1.f);
      }

      ++finished;
      if (finished % 20 == 0) {
        printProgress((float)finished / (height * width));
      }
    }
#endif
  printProgress(1);
}