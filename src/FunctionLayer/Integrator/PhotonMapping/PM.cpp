#include "PM.h"
#include "nanoflann/nanoflann.hpp"
#include <FunctionLayer/Material/Material.h>
#include <mutex>
#include <tbb/tbb.h>
struct Photon {
  Spectrum weight;
  Point3f position;
  Vector3f wi;
  std::shared_ptr<Phase> phase;
  std::shared_ptr<BSDF> bsdf;
};

class PhotonMap {
  std::mutex mtx;
  std::vector<Photon> photons;

public:
  PhotonMap(size_t size) { photons.reserve(size); }

  void emplace_back(Photon p) {
    std::lock_guard<std::mutex> lk_gd(mtx);
    photons.emplace_back(p);
  }

  size_t kdtree_get_point_count() const { return photons.size(); }

  float kdtree_get_pt(const size_t idx, const size_t dim) const {
    return photons[idx].position[dim];
  }

  // Just return false
  template <class BBox> bool kdtree_get_bbox(BBox & /*bb*/) const {
    return false;
  }

  const Photon &get_photon(size_t index) const { return photons[index]; }
};

using my_kd_tree_t = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<float, PhotonMap>, PhotonMap, 3 /* dim */
    >;

PhotonMapping::PhotonMapping(const Json &json) : Integrator(json) {
  maxDepth = fetchRequired<int>(json, "maxPathLength");
  k = fetchRequired<int>(json, "k");
  N = fetchRequired<int>(json, "N");
}

void PhotonMapping::render(const Camera &camera, const Scene &scene,
                           std::shared_ptr<Sampler> sampler, int spp) const {
  //* Generate N photons and store them in kd-tree
  PhotonMap photonMap(maxDepth * N);
  std::cout << "Generate photons ...\n";
  std::atomic<int> n = 0;
  tbb::parallel_for(tbb::blocked_range<size_t>(0, N), [&](const tbb::
                                                              blocked_range<
                                                                  size_t> &r) {
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

        //*------ PhotonRay Random Walk ------
        for (int depth = 0; depth < maxDepth; ++depth) {
          auto si = scene.rayIntersect(photonRay);

          const Medium *medium = photonRay.medium; // TODO

          if (!si /* Terminate */)
            break;

          auto material = si->shape->material;
          auto bsdf = material->computeBSDF(*si);

          if (!bsdf /* Skip the empty surface */) {
            photonRay = Ray{si->position, photonRay.direction, 1e-4f, FLT_MAX};
            --depth;
            continue;
          }

          if (bsdf->type ==
              BSDFType::Specular /* Spwan at specular surface */) {
            BSDFSampleResult result =
                bsdf->sample(-photonRay.direction, sampler->next2D());
            beta *= result.weight;
            photonRay = Ray{si->position, result.wi, 1e-4f, FLT_MAX};
            if (beta.isZero())
              break;
            continue;
          }

          // Store the photon on diffuse surface
          Photon photon{beta, si->position, -photonRay.direction, nullptr,
                        bsdf};
          photonMap.emplace_back(photon);

          // Spwan the photon ray
          BSDFSampleResult result =
              bsdf->sample(-photonRay.direction, sampler->next2D());
          beta *= result.weight;
          if (beta.isZero())
            break;
          photonRay = Ray{si->position, result.wi, 1e-4f, FLT_MAX};
        }
      }
      if (++n % 100 == 0) {
        printProgress((float)n / N);
      }
    }
  });
  printProgress(1.f);
  std::cout << "\n";

  my_kd_tree_t index(3, photonMap, 128);

  int width = camera.film->size[0];
  int height = camera.film->size[1];

  n = 0;
  std::cout << "Gathering photons ...\n";
  for (int i = 0; i < spp; ++i) {
    tbb::parallel_for(
        tbb::blocked_range2d<size_t>(0, width, 0, height),
        [&](const tbb::blocked_range2d<size_t> &r) {
          for (int row = r.rows().begin(); row != r.rows().end(); ++row) {
            for (int col = r.cols().begin(); col != r.cols().end(); ++col) {
              Vector2f NDC{(float)row / width, (float)col / height};
              Ray ray = camera.sampleRayDifferentials(
                  CameraSample{sampler->next2D()}, NDC);

              Spectrum Li(.0f);
              Spectrum beta(1.f);
              bool specularBounce = false;

              for (int depth = 0; depth < maxDepth; ++depth) {
                auto si = scene.rayIntersect(ray);

                // TODO add medium
                if (!si /* Terminate */)
                  break;

                auto material = si->shape->material;
                auto bsdf = material->computeBSDF(*si);

                if (!bsdf /* Skip empty surface */) {
                  ray = Ray{si->position, ray.direction, 1e-4f, FLT_MAX};
                  --depth;
                  continue;
                }

                if (specularBounce || depth == 0) {
                  if (auto light = si->shape->light; light) {
                    Li += beta * light->evaluateEmission(*si, -ray.direction);
                  }
                }

                if (bsdf->type ==
                    BSDFType::Specular /* Spwan specular bounce */) {
                  specularBounce = true;

                  BSDFSampleResult result =
                      bsdf->sample(-ray.direction, sampler->next2D());
                  beta *= result.weight;
                  ray = Ray{si->position, result.wi, 1e-4f, FLT_MAX};
                  continue;
                }

                // Search the k-nearest photon for diffuse hitpoint

                float query_point[3] = {si->position[0], si->position[1],
                                        si->position[2]};
                std::vector<uint32_t> ret_index(k);
                std::vector<float> out_dist_sqr(k);

                size_t num = index.knnSearch(query_point, k, &ret_index[0],
                                             &out_dist_sqr[0]);

                Spectrum phi(.0f);
                float r_max = .0f;
                for (int j = 0; j < num; ++j) {
                  auto photon = photonMap.get_photon(ret_index[j]);
                  auto f = photon.bsdf->f(-ray.direction, photon.wi);
                  phi += beta * photon.weight * f;
                  r_max = std::max(r_max, out_dist_sqr[j]);
                }

                if (!phi.isZero())
                  Li += phi / ((float)N * PI * r_max);

                if (++n % 100 == 0) {
                  printProgress((float)n / (width * height * spp));
                }

                camera.film->deposit({row, col}, Li, 1.f);
                break;
              }
            }
          }
        });
  }
  printProgress(1.f);
}

REGISTER_CLASS(PhotonMapping, "pm")