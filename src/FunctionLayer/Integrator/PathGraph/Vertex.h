#pragma once
#include <CoreLayer/ColorSpace/Spectrum.h>
#include <CoreLayer/Math/Math.h>
#include <FunctionLayer/Film/Film.h>
#include <FunctionLayer/Light/Light.h>
#include <FunctionLayer/Material/BxDF/BSDF.h>
#include <FunctionLayer/Shape/Shape.h>
#include <memory>
#include <mutex>

enum class VertexType { SURFACE, CAMERA, LIGHT };

struct Vertex;

// Always from light towards camera
struct DirectionalEdge {
  DirectionalEdge(std::shared_ptr<Vertex> f, std::shared_ptr<Vertex> t,
                  Spectrum w)
      : from(f.get()), to(t.get()), weight(w), radiance(Spectrum(.0f)) {}

  const Vertex *from = nullptr, *to = nullptr;
  Spectrum weight;
  Spectrum radiance;
};

struct Vertex {
  Vertex(Point3f p, Vector3f n, VertexType t)
      : position(p), normal(n), type(t) {}

  VertexType type;
  Point3f position;
  Vector3f normal;
  std::shared_ptr<DirectionalEdge> next;
  std::vector<std::shared_ptr<Vertex>> prevVertices;
  bool visited = false; // For tarverse
};

struct SurfaceVertex : public Vertex {
  SurfaceVertex(Point3f p, Vector3f n) : Vertex(p, n, VertexType::SURFACE) {}

  const Shape *shape = nullptr;
  std::shared_ptr<BSDF> bsdf;
  Vector3f wo;
};

struct CameraVertex : public Vertex {
  CameraVertex() : Vertex(Point3f(), Vector3f(), VertexType::CAMERA) {}
};

struct LightVertex : public Vertex {
  LightVertex(Point3f p, Vector3f n) : Vertex(p, n, VertexType::LIGHT) {}

  const Shape *shape = nullptr;
  const Light *light = nullptr;
  Spectrum energy;
};

struct LightPath {
  LightPath(Vector2i _pixelLoc)
      : pixelLoc(_pixelLoc), root(std::make_shared<CameraVertex>()), cur(root) {
  }

  void attachInfLightVertex(Spectrum weight,
                            std::shared_ptr<InfiniteLight> light,
                            const Ray &ray);

  void attachLightVertex(Spectrum weight, std::shared_ptr<Light> light,
                         const Intersection &its, const Ray &ray);

  void attachLightVertex(Spectrum weight, std::shared_ptr<Light> light,
                         const LightSampleResult &result, const Ray &ray);

  void attachSurfaceVertex(Spectrum weight, const Intersection &its,
                           std::shared_ptr<BSDF> bsdf, Vector3f wo);

  Spectrum gatherRadiance() const;

public:
  Vector2i pixelLoc;

protected:
  std::shared_ptr<Vertex>
  currentVertex(std::shared_ptr<SurfaceVertex> surfaceVertex = nullptr);

  std::shared_ptr<CameraVertex> root;

  std::shared_ptr<Vertex> cur;
};

struct PathGraph {
  PathGraph(Vector2i resolution, int _spp) : spp(_spp) {
    paths.reserve(resolution[0] * resolution[1] * spp);
  }

  void addPath(LightPath path) {
    std::lock_guard<std::mutex> lk(mtx);
    paths.emplace_back(path);
  }

  void toFilm(std::shared_ptr<Film> film) const;

public:
  std::vector<LightPath> paths;
  std::mutex mtx;
  int spp;
};