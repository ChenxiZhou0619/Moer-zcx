#include <CoreLayer/Math/Math.h>
#include <FunctionLayer/Camera/Pinhole.h>
#include <FunctionLayer/Integrator/Integrator.h>
#include <FunctionLayer/Sampler/Sampler.h>
#include <FunctionLayer/Scene/Scene.h>
#include <FunctionLayer/Texture/Mipmap.h>
#include <ResourceLayer/Factory.h>
#include <ResourceLayer/FileUtil.h>
#include <ResourceLayer/Image.h>
#include <ResourceLayer/JsonUtil.h>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <regex>
#include <stdio.h>
#include <tbb/tbb.h>

int main(int argc, char **argv) {
  const std::string sceneDir = std::string(argv[1]);
  FileUtil::setWorkingDirectory(sceneDir);
  std::string sceneJsonPath = FileUtil::getFullPath("scene.json");
  std::ifstream fstm(sceneJsonPath);
  Json json = Json::parse(fstm);
  auto camera = Factory::construct_class<Camera>(json["camera"]);
  auto scene = std::make_shared<Scene>(json["scene"]);
  auto integrator = Factory::construct_class<Integrator>(json["integrator"]);
  if (scene->sceneMedium) {
    camera->medium = scene->sceneMedium.get();
    integrator->medium = scene->sceneMedium.get();
  }
  auto sampler = Factory::construct_class<Sampler>(json["sampler"]);
  int spp = sampler->xSamples * sampler->ySamples;

  auto start = std::chrono::system_clock::now();
  integrator->render(*camera, *scene, sampler, spp);
  auto end = std::chrono::system_clock::now();

  printf("\nRendering costs %.2fs\n",
         (std::chrono::duration_cast<std::chrono::milliseconds>(end - start))
                 .count() /
             1000.f);

  std::string outputDir = FileUtil::getFullPath("result");
  if (!std::filesystem::exists(outputDir)) {
    std::filesystem::create_directory(outputDir);
  }

  std::string outputName =
      outputDir + "/" + fetchRequired<std::string>(json["output"], "filename");

  auto path = saveImage(outputName.c_str(), camera->film->image, false);
  std::cout << "Save at " << path << std::endl;
}
