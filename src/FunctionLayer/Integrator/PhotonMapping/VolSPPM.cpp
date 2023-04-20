#include "VolSPPM.h"
#include <FunctionLayer/Acceleration/AABB.h>
#include <FunctionLayer/Material/Material.h>
#include <atomic>
#include <tbb/tbb.h>
void VolSPPM::setRayMedium(Vector3f direction, Vector3f normal,
                           std::shared_ptr<Material> material, Ray *ray) const {
  bool towardsInner = dot(direction, normal) < .0f;
  ray->medium = towardsInner ? material->getMedium() : nullptr;
}

//* For each pixel each iteration, we store the its corresponding information
struct VolSPPMPixel {
  Spectrum Phi;           // Flux of current iteration
  std::atomic<int> M = 0; // Photons collected by this vp this iteration
  float N = 0;            // Photons collected totally
  Spectrum tau{.0f};      // Flux estimated by previous iterations
  Spectrum Ld{.0f};       // Le and Ld which is not estimated by sppm
  float radius;           // The search radius of current iteration

  mutable std::mutex mtx;

  void addPhi(Spectrum phi) {
    std::lock_guard<std::mutex> lg(mtx);
    Phi += phi;
  }

  struct VisiblePoint {
    Point3f position;             // Position of the vp
    Vector3f normal;              // Normal at the vp
    std::shared_ptr<BSDF> bsdf;   // BSDF at vp
    std::shared_ptr<Phase> phase; // Phase at vp
    Spectrum alpha{.0f};          // Weight of vp to  pixel
    Vector3f wo;                  // Omega_o
    bool isSurface;               // Whether the vp is on surface
  } vp;
};

struct VolSPPMPixelListNode {
  VolSPPMPixelListNode *next = nullptr;
  VolSPPMPixel *pixel = nullptr;

  ~VolSPPMPixelListNode() {
    if (next)
      delete next;
  }
};

VolSPPM::VolSPPM(const Json &json) : Integrator(json) {
  maxDepth = fetchOptional(json, "maxPathLength", 5);
  photonsPerIteration = fetchRequired<int>(json, "photonsPerIteration");
  searchRadius = fetchOptional(json, "searchRadius", 0.005f);
}

bool VolToGrid(Point3f position, const AABB &bounds, int gridRes[3],
               int gridCoord[3]) {

  bool inGrid = true;

  Point3f uCoord = bounds.UniformCoord(position);
  for (int i = 0; i < 3; ++i) {
    int iCoord = (int)(uCoord[i] * gridRes[i]);
    if (iCoord < 0 || iCoord >= gridRes[i])
      inGrid = false;
    gridCoord[i] = clamp<int>(iCoord, 0, gridRes[i] - 1);
  }
  return inGrid;
}

unsigned int Volhash(int x, int y, int z, int hashSize) {
  return (unsigned int)((x * 73856093) ^ (y * 19349663) ^ (z * 83492791)) %
         hashSize;
}

float VolDistanceSquare(Point3f p1, Point3f p2) {
  float distance = (p1 - p2).length();
  return distance * distance;
}

void VolSPPM::render(const Camera &camera, const Scene &scene,
                     std::shared_ptr<Sampler> sampler, int spp) const {
  // 1. Sample camera ray, store the visible points in a hash grid
  // 2. Generate photons, for each photon, find its neighbor vps, add the
  // contribution to them
  // 3. Continue iterations

  // // TODO Sample visible points
  // // TODO Spatial hash grid
  // // TODO Generate photons and corresponding beta
  // // TODO Gather contribution
  // // TODO Iterations
  // // TODO Finally, reconstruct image

  //! Notice SPPM doesn't support environment map now

  // Construct SPPMPixel for each pixels
  int width = camera.film->size[0];
  int height = camera.film->size[1];
  int nPixels = width * height;
  std::unique_ptr<VolSPPMPixel[]> pixels(new VolSPPMPixel[nPixels]);

  for (int i = 0; i < nPixels; ++i) {
    pixels[i].radius = searchRadius;
  }

  int sum = 0;
  int iterations = spp;

  for (int itr = 0; itr < iterations; ++itr) {
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
              VolSPPMPixel *pixel = &pixels[row + col * width];
              bool specular_bounce = false;

              //* ----- Random walk -----
              for (int depth = 0; depth < maxDepth; ++depth) {
                auto si = scene.rayIntersect(ray);
                Vector3f wo = -ray.direction;

                //* Sample medium if exists
                const Medium *medium = ray.medium;
                MediumIntersection mi;
                bool mediumInteraction = false;

                if (medium) {
                  Spectrum Tr;
                  float pdf, tmax = si ? si->distance : FLT_MAX;
                  mediumInteraction = medium->Sample_RegularTracking(
                      ray, tmax, sampler->next2D(), &mi, &Tr, &pdf);
                  alpha *= Tr / pdf;
                }
                if (mediumInteraction) {
                  alpha *= mi.mp.sigma_s;
                  pixel->vp = {mi.position, Vector3f(), nullptr, mi.mp.phase,
                               alpha,       wo,         false};
                  break;
                }

                if (!si /* Terminate if escape the scene */) {
                  //* Just break cuz no support for environmentlight
                  break;
                }

                //* Account for emission term
                if (depth == 0 || specular_bounce) {
                  if (auto light = si->shape->light; light) {
                    pixel->Ld += alpha * light->evaluateEmission(*si, wo);
                  }
                }

                auto material = si->shape->material;
                auto bsdf = material->computeBSDF(*si);

                //* Skip empty surface
                if (!bsdf) {
                  ray = Ray{si->position, ray.direction, 1e-4f, FLT_MAX};
                  setRayMedium(ray.direction, si->normal, material, &ray);
                  --depth;
                  continue;
                }

                //* If hit the diffuse surface, terminate and record vp
                if (bsdf->type == BSDFType::Diffuse) {
                  pixel->vp = {si->position, si->normal, bsdf, nullptr,
                               alpha,        wo,         true}; // TODO
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

    const int hashSize = nPixels;
    std::vector<std::atomic<VolSPPMPixelListNode *>> grid(hashSize);
    AABB vp_bound;
    int gridRes[3];

    //* Organize visible points in hash grid
    {
      float max_radius = .0f;
      for (int i = 0; i < nPixels; ++i) {
        const auto &pixel = pixels[i];
        if (pixel.vp.alpha.isZero())
          continue;
        max_radius = std::max(max_radius, pixel.radius);
        AABB search_box{pixel.vp.position - Vector3f{pixel.radius},
                        pixel.vp.position + Vector3f{pixel.radius}};
        vp_bound = vp_bound.Union(search_box);
      }

      Vector3f diag = vp_bound.pMax - vp_bound.pMin;
      float diagMax = MaxComponent(diag);
      int gridResBase = (int)(diagMax / max_radius);

      for (int i = 0; i < 3; ++i) {
        gridRes[i] = std::max((int)(gridResBase * diag[i] / diagMax), 1);
      }

      tbb::parallel_for(tbb::blocked_range<size_t>(0, nPixels),
                        [&](const tbb::blocked_range<size_t> &r) {
                          for (int i = r.begin(); i != r.end(); ++i) {
                            VolSPPMPixel *pixel = &pixels[i];
                            if (pixel->vp.alpha.isZero())
                              continue;
                            float radius = pixel->radius;

                            int pMin[3], pMax[3];
                            VolToGrid(pixel->vp.position - Vector3f{radius},
                                      vp_bound, gridRes, pMin);
                            VolToGrid(pixel->vp.position + Vector3f{radius},
                                      vp_bound, gridRes, pMax);

                            for (int x = pMin[0]; x <= pMax[0]; ++x)
                              for (int y = pMin[1]; y <= pMax[1]; ++y)
                                for (int z = pMin[2]; z <= pMax[2]; ++z) {
                                  int h = Volhash(x, y, z, hashSize);
                                  VolSPPMPixelListNode *node =
                                      new VolSPPMPixelListNode();
                                  node->pixel = pixel;
                                  node->next = grid[h];
                                  while (grid[h].compare_exchange_weak(
                                             node->next, node) == false)
                                    ;
                                }
                          }
                        });
    }

    //* Trace photons and compute contribution
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, photonsPerIteration),
        [&](const tbb::blocked_range<size_t> &r) {
          for (int i = r.begin(); i != r.end(); ++i) {
            float pdfLight;
            auto light = scene.sampleLight(sampler->next1D(), &pdfLight);
            if (light && pdfLight != .0f) {
              float pdfPhotonRay;
              Ray photonRay;
              Spectrum Le;
              Vector3f nLight;
              light->sampleLe(sampler->next2D(), sampler->next2D(), &photonRay,
                              &pdfPhotonRay, &Le, &nLight);
              Spectrum beta = std::abs(dot(photonRay.direction, nLight)) * Le /
                              (pdfPhotonRay * pdfLight);
              if (beta.isZero())
                continue;

              //* ------------------------------
              //* ----- Photon Random Walk -----
              //* ------------------------------

              for (int depth = 0; depth < maxDepth; ++depth) {
                auto si = scene.rayIntersect(photonRay);

                //* Sample medium if exists
                const Medium *medium = photonRay.medium;
                MediumIntersection mi;
                bool mediumInteraction = false;

                if (medium) {
                  Spectrum Tr;
                  float pdf, tmax = si ? si->distance : FLT_MAX;
                  mediumInteraction = medium->Sample_RegularTracking(
                      photonRay, tmax, sampler->next2D(), &mi, &Tr, &pdf);
                  // TODO Beta change
                }

                if (mediumInteraction) {
                  int gridCoord[3];
                  if (VolToGrid(mi.position, vp_bound, gridRes, gridCoord)) {
                    int h = Volhash(gridCoord[0], gridCoord[1], gridCoord[2],
                                    hashSize);
                    for (auto *node = grid[h].load(std::memory_order_relaxed);
                         node != nullptr; node = node->next) {
                      auto pixel = node->pixel;
                      float radius = pixel->radius;
                      if ((VolDistanceSquare(pixel->vp.position, mi.position) <
                           radius * radius) &&
                          !pixel->vp.isSurface) {
                        Spectrum f = pixel->vp.phase->f(pixel->vp.wo,
                                                        -photonRay.direction);
                        Spectrum phi = beta * f;
                        pixel->addPhi(phi);
                        ++(pixel->M);
                      }
                    }
                  }

                  //* Sample phase to spwan the ray
                  PhaseSampleResult result = mi.mp.phase->sample(
                      -photonRay.direction, sampler->next2D());

                  beta *= result.weight;
                  if (beta.isZero())
                    break;
                  photonRay = Ray{si->position, result.wi, 1e-4f, FLT_MAX};
                  continue;
                }

                if (!si)
                  break;

                {
                  auto material = si->shape->material;
                  auto bsdf = material->computeBSDF(*si);
                  if (!bsdf /* Skip the empty surface */) {
                    photonRay =
                        Ray{si->position, photonRay.direction, 1e-4f, FLT_MAX};
                    setRayMedium(photonRay.direction, si->normal, material,
                                 &photonRay);
                    --depth;
                    continue;
                  }
                }

                int gridCoord[3];
                if (VolToGrid(si->position, vp_bound, gridRes, gridCoord)) {
                  int h = Volhash(gridCoord[0], gridCoord[1], gridCoord[2],
                                  hashSize);
                  for (auto *node = grid[h].load(std::memory_order_relaxed);
                       node != nullptr; node = node->next) {
                    auto pixel = node->pixel;
                    float radius = pixel->radius;
                    if ((VolDistanceSquare(pixel->vp.position, si->position) <
                         radius * radius) &&
                        pixel->vp.isSurface) {
                      Spectrum f =
                          pixel->vp.bsdf->f(pixel->vp.wo, -photonRay.direction);
                      Spectrum phi = beta * f;
                      pixel->addPhi(phi);
                      ++(pixel->M);
                    }
                  }
                }

                //* Spwan the photon ray
                auto material = si->shape->material;
                auto bsdf = material->computeBSDF(*si);

                BSDFSampleResult result =
                    bsdf->sample(-photonRay.direction, sampler->next2D());
                beta *= result.weight;
                if (beta.isZero())
                  break;
                photonRay = Ray{si->position, result.wi, 1e-4f, FLT_MAX};
                setRayMedium(photonRay.direction, si->normal, material,
                             &photonRay);
              }
            }
            if ((++sum % 100) == 0) {
              printProgress((float)sum / (photonsPerIteration * iterations));
            }
          }
        });

    for (int i = 0; i < nPixels; ++i) {
      if (grid[i]) {
        delete grid[i];
      }
    }

    //* Update pixel information
    tbb::parallel_for(
        tbb::blocked_range<size_t>(0, nPixels),
        [&](const tbb::blocked_range<size_t> &r) {
          for (int i = r.begin(); i != r.end(); ++i) {
            VolSPPMPixel *pixel = &pixels[i];
            if (pixel->M > 0) {
              constexpr float gamma = 2.f / 3.f;
              float Nnew = pixel->N + gamma * pixel->M;
              float Rnew =
                  pixel->radius * fm::sqrt(Nnew / (pixel->N + pixel->M));
              pixel->tau = pixel->tau + pixel->vp.alpha * pixel->Phi;
              pixel->tau *= (Rnew * Rnew) / (pixel->radius * pixel->radius);
              pixel->N = Nnew;
              pixel->radius = Rnew;
              pixel->M = 0;
              pixel->Phi = Spectrum(.0f);
            }
            pixel->vp.alpha = .0f;
            pixel->vp.bsdf = nullptr;
            pixel->vp.phase = nullptr;
          }
        });
  }
  printProgress(1.f);

  //* Reconstruct the image
  tbb::parallel_for(tbb::blocked_range<size_t>(0, nPixels),
                    [&](const tbb::blocked_range<size_t> &r) {
                      for (int i = r.begin(); i != r.end(); ++i) {
                        int x = i % width;
                        int y = i / width;
                        Spectrum L = pixels[i].Ld / iterations;
                        float radius = pixels[i].radius;
                        L += pixels[i].tau / (iterations * photonsPerIteration *
                                              PI * radius * radius);
                        camera.film->deposit({x, y}, L, 1.f);
                      }
                    });
}

REGISTER_CLASS(VolSPPM, "volsppm")