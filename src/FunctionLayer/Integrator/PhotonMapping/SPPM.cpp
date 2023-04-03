#include "SPPM.h"
#include <FunctionLayer/Acceleration/AABB.h>
#include <FunctionLayer/Material/Material.h>
#include <atomic>
#include <tbb/tbb.h>
//* For each pixel each iteration, we store the its corresponding information
struct SPPMPixel {
  std::atomic<Spectrum> phi{.0f}; // Flux of current iteration
  std::atomic<int> M = 0;         // Photons collected by this vp this iteration
  Spectrum tau{.0f};              // Flux estimated by previous iterations
  Spectrum Ld{.0f};               // Le and Ld which is not estimated by sppm
  float radius;                   // The search radius of current iteration

  struct VisiblePoint {
    Point3f position;           // Position of the vp
    Vector3f normal;            // Normal at the vp
    std::shared_ptr<BSDF> bsdf; // BSDF at vp
    Spectrum alpha{.0f};        // Weight of vp to  pixel
  } vp;
};

SPPM::SPPM(const Json &json) : Integrator(json) {
  //
}

void SPPM::render(const Camera &camera, const Scene &scene,
                  std::shared_ptr<Sampler> sampler, int spp) const {
  // 1. Sample camera ray, store the visible points in a hash grid
  // 2. Generate photons, for each photon, find its neighbor vps, add the
  // contribution to them
  // 3. Continue iterations

  // TODO Sample visible points
  // TODO Spatial hash grid
  // TODO Generate photons and corresponding beta
  // TODO Gather contribution
  // TODO Iterations
  // TODO Finally, reconstruct image

  //! Notice SPPM doesn't support environment map now

  // Construct SPPMPixel for each pixels
  int width = camera.film->size[0];
  int height = camera.film->size[1];
  int nPixels = width * height;
  std::unique_ptr<SPPMPixel[]> pixels(new SPPMPixel[nPixels]);

  tbb::parallel_for(
      tbb::blocked_range2d<size_t>(0, width, 0, height),
      [&](const tbb::blocked_range2d<size_t> &r) {
        for (int row = r.rows().begin(); row != r.rows().end(); ++row)
          for (int col = r.cols().begin(); col != r.cols().end(); ++col) {
            //* ------------------------
            //* 1. Sample visible points
            //* ------------------------

            Vector2f NDC{(float)row / width, (float)col / height};
            Ray ray = camera.sampleRayDifferentials(
                CameraSample{sampler->next2D()}, NDC);
            Spectrum alpha(1.f);
            SPPMPixel *pixel = &pixels[row + col * width];
            bool specular_bounce = false;

            //* ----- Random walk -----
            for (int depth = 0; depth < maxDepth; ++depth) {
              auto si = scene.rayIntersect(ray);
              Vector3f wo = -ray.direction;
              if (!si /* Terminate if escape the scene */) {
                //* Just break cuz no support for environmentlight
                break;
              }

              auto material = si->shape->material;
              auto bsdf = material->computeBSDF(*si);

              //* Account for emission term
              if (depth == 0 || specular_bounce) {
                if (auto light = si->shape->light; light) {
                  pixel->Ld += alpha * light->evaluateEmission(*si, wo);
                }
              }

              //* If hit the diffuse surface, terminate and record vp
              if (bsdf->type == BSDFType::Diffuse) {
                pixel->vp = {si->position, si->normal, bsdf, alpha};
                break;
              } else /* Hit specular surface */ {
                BSDFSampleResult result = bsdf->sample(wo, sampler->next2D());
                alpha *= result.weight;
                if (alpha.isZero())
                  break;
                ray = Ray{si->position, result.wi, 1e-4f, FLT_MAX};
                specular_bounce = true;
              }
            }
          }
      });

  //* Organize visible points in hash grid
  {
    float max_radius = .0f;
    AABB vp_bound;
    for (int i = 0; i < nPixels; ++i) {
      const auto &pixel = pixels[i];
      if (pixel.vp.alpha.isZero())
        continue;
      max_radius = std::max(max_radius, pixel.radius);
      AABB search_box{pixel.vp.position - Vector3f{pixel.radius},
                      pixel.vp.position + Vector3f{pixel.radius}};
      vp_bound.Union(search_box);
    }

    int gridRes[3];
    Vector3f diag = vp_bound.pMax - vp_bound.pMin;
    float diagMax = MaxComponent(diag);
    int gridResBase = (int)(diagMax / max_radius);

    for (int i = 0; i < 3; ++i) {
      gridRes[i] = std::max((int)(gridResBase * diag[i] / diagMax), 1);
    }
  }
  //* Trace photons and compute contribution

  //* Reconstruct the image
  tbb::parallel_for(tbb::blocked_range<size_t>(0, nPixels),
                    [&](const tbb::blocked_range<size_t> &r) {
                      for (int i = r.begin(); i != r.end(); ++i) {
                        int x = i % width;
                        int y = i / width;
                        Spectrum L; // TODO compute L
                        camera.film->deposit({x, y}, L, 1.f);
                      }
                    });
}

REGISTER_CLASS(SPPM, "sppm")