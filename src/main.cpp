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
#include <fstream>
#include <openvdb/openvdb.h>
#include <regex>
#include <stdio.h>
#include <tbb/tbb.h>

#include <nanovdb/NanoVDB.h>
#include <nanovdb/util/IO.h>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

inline void printProgress(float percentage) {
  int val = (int)(percentage * 100);
  int lpad = (int)(percentage * PBWIDTH);
  int rpad = PBWIDTH - lpad;
  printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
  fflush(stdout);
}

int main(int argc, char **argv) {
  const std::string sceneDir = std::string(argv[1]);
  FileUtil::setWorkingDirectory(sceneDir);
  std::string sceneJsonPath = FileUtil::getFullPath("scene.json");
  std::ifstream fstm(sceneJsonPath);
  Json json = Json::parse(fstm);
  auto camera = Factory::construct_class<Camera>(json["camera"]);
  auto scene = std::make_shared<Scene>(json["scene"]);
  auto integrator = Factory::construct_class<Integrator>(json["integrator"]);
  auto sampler = Factory::construct_class<Sampler>(json["sampler"]);
  int spp = sampler->xSamples * sampler->ySamples;
  int width = camera->film->size[0], height = camera->film->size[1];

  auto start = std::chrono::system_clock::now();

  int finished = 0;
  tbb::parallel_for(
      tbb::blocked_range2d<size_t>(0, width, 0, height),
      [&](const tbb::blocked_range2d<size_t> &r) {
        for (int row = r.rows().begin(); row != r.rows().end(); ++row)
          for (int col = r.cols().begin(); col != r.cols().end(); ++col) {
            Vector2f NDC{(float)row / width, (float)col / height};
            Spectrum li(.0f);
            for (int i = 0; i < spp; ++i) {
              Ray ray = camera->sampleRayDifferentials(
                  CameraSample{sampler->next2D()}, NDC);
              li += integrator->li(ray, *scene, sampler);
            }
            camera->film->deposit({row, col}, li / spp);

            ++finished;
            if (finished % 20 == 0) {
              printProgress((float)finished / (height * width));
            }
          }
      });
  printProgress(1);

  auto end = std::chrono::system_clock::now();

  printf("\nRendering costs %.2fs\n",
         (std::chrono::duration_cast<std::chrono::milliseconds>(end - start))
                 .count() /
             1000.f);

  //* 目前支持输出为png/hdr两种格式
  std::string outputName =
      fetchRequired<std::string>(json["output"], "filename");
  if (std::regex_match(outputName, std::regex("(.*)(\\.png)"))) {
    camera->film->savePNG(outputName.c_str());
  } else if (std::regex_match(outputName, std::regex("(.*)(\\.hdr)"))) {
    camera->film->saveHDR(outputName.c_str());
  } else {
    std::cout << "Only support output as PNG/HDR\n";
  }
}
