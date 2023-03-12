#include "Vertex.h"
#include <FunctionLayer/Light/Light.h>
#include <stack>
void LightPath::attachInfLightVertex(Spectrum weight,
                                     std::shared_ptr<InfiniteLight> light,
                                     const Ray &ray) {
  auto cur = currentVertex();
  // Construct New Vertex
  std::shared_ptr<LightVertex> lv =
      std::make_shared<LightVertex>(Point3f(.0f), Vector3f(.0f));
  lv->light = light.get();
  lv->shape = nullptr;
  lv->energy = light->evaluateEmission(ray);
  // Construct New Edge
  std::shared_ptr<DirectionalEdge> toCur =
      std::make_shared<DirectionalEdge>(lv, cur, weight);
  lv->next = toCur;
  cur->prevVertices.push_back(lv);
}

void LightPath::attachLightVertex(Spectrum weight, std::shared_ptr<Light> light,
                                  const Intersection &its, const Ray &ray) {
  auto cur = currentVertex();
  // Construct New Vertex
  std::shared_ptr<LightVertex> lv =
      std::make_shared<LightVertex>(its.position, its.normal);
  lv->light = light.get();
  lv->shape = its.shape;
  lv->energy = light->evaluateEmission(its, -ray.direction);
  // Construct New Edge
  std::shared_ptr<DirectionalEdge> toCur =
      std::make_shared<DirectionalEdge>(lv, cur, weight);
  lv->next = toCur;
  cur->prevVertices.push_back(lv);
}

void LightPath::attachLightVertex(Spectrum weight, std::shared_ptr<Light> light,
                                  const LightSampleResult &result,
                                  const Ray &ray) {
  auto cur = currentVertex();
  // Construct New Vertex
  std::shared_ptr<LightVertex> lv =
      std::make_shared<LightVertex>(ray.at(result.distance), result.normal);
  lv->light = light.get();
  lv->shape = nullptr; // TODO delete this?
  lv->energy = result.energy;
  // Construct New Edge
  std::shared_ptr<DirectionalEdge> toCur =
      std::make_shared<DirectionalEdge>(lv, cur, weight);
  cur->prevVertices.push_back(lv);
}

void LightPath::attachSurfaceVertex(Spectrum weight, const Intersection &its,
                                    std::shared_ptr<BSDF> bsdf, Vector3f wo) {
  auto cur = currentVertex();
  // Construct New Vertex
  std::shared_ptr<SurfaceVertex> sv =
      std::make_shared<SurfaceVertex>(its.position, its.normal);
  sv->shape = its.shape;
  sv->bsdf = bsdf;
  sv->wo = wo;
  // Construct New Edge
  std::shared_ptr<DirectionalEdge> toCur =
      std::make_shared<DirectionalEdge>(sv, cur, weight);
  sv->next = toCur;
  cur->prevVertices.push_back(sv);
  currentVertex(sv);
}

Spectrum LightPath::gatherRadiance() const {
  std::stack<std::shared_ptr<Vertex>> stack;
  std::shared_ptr<Vertex> node = root;
  stack.push(node);

  while (!stack.empty()) {
    node = stack.top();
    stack.pop();

    if (node->prevVertices.empty()) {
      // End point
      if (node->type == VertexType::LIGHT) {
        auto tmp = std::static_pointer_cast<LightVertex>(node);
        node->next->radiance = tmp->energy;
      }
      node->visited = true;
      continue;
    }

    if (!node->visited && !node->prevVertices.empty()) {
      node->visited = true;
      stack.push(node);
      for (auto prev : node->prevVertices) {
        stack.push(prev);
      }
      continue;
    }

    if (node->visited) {
      Spectrum localSum(.0f);
      for (auto prev : node->prevVertices) {
        localSum += prev->next->radiance * prev->next->weight;
      }
      if (node->type != VertexType::CAMERA)
        node->next->radiance = localSum;
    }
  }

  Spectrum spectrum(.0f);
  for (auto prev : root->prevVertices) {
    spectrum += prev->next->radiance * prev->next->weight;
  }
  return spectrum;
}

std::shared_ptr<Vertex>
LightPath::currentVertex(std::shared_ptr<SurfaceVertex> surfaceVertex) {
  if (surfaceVertex)
    cur = surfaceVertex;
  return cur;
}

void PathGraph::toFilm(std::shared_ptr<Film> film) const {
  for (const auto &path : paths) {
    //
  }
}