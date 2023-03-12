#include "Integrator.h"

#include <tbb/tbb.h>

void PixelIntegrator::render(const Camera &camera, const Scene &scene,
                             std::shared_ptr<Sampler> sampler, int spp) const {
  int width = camera.film->size[0], height = camera.film->size[1], finished = 0;
  tbb::parallel_for(
      tbb::blocked_range2d<size_t>(0, width, 0, height),
      [&](const tbb::blocked_range2d<size_t> &r) {
        for (int row = r.rows().begin(); row != r.rows().end(); ++row)
          for (int col = r.cols().begin(); col != r.cols().end(); ++col) {
            Vector2f NDC{(float)row / width, (float)col / height};
            Spectrum spectrum(.0f);
            for (int i = 0; i < spp; ++i) {
              Ray ray = camera.sampleRayDifferentials(
                  CameraSample{sampler->next2D()}, NDC);
              spectrum += li(ray, scene, sampler);
            }
            camera.film->deposit({row, col}, spectrum / spp);

            ++finished;
            if (finished % 20 == 0) {
              printProgress((float)finished / (height * width));
            }
          }
      });
  printProgress(1);
}