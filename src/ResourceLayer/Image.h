#pragma once
#include <CoreLayer/Math/Math.h>
#include <memory>

// TODO 实现一个模板类，提供图片多通道、多数据类型的支持
//* 目前只支持三通道（即RGB）的图片格式

//* Image的数据组织格式如下
//* [0, 0] ----X-----
//* |               |
//* |               |
//* Y        [x, y] |
//* -----------------
//* 坐标[0, 0]对应图片左上角
//* 坐标[x, y]对应图片第x列，第y行

class Image {
public:
  Image() = delete;

  Image(Vector2i _size) : size(_size) {
    data = new float[_size[0] * _size[1] * channels]();
    weight = new float[_size[0] * _size[1]]();
  }

  Image(Vector2i _size, float *_data) : size(_size), data(_data) {
    weight = new float[_size[0] * _size[1]]();
    for (int i = 0; i < size[0] * size[1]; ++i)
      weight[i] = 1.f;
  }

  ~Image() {
    delete[] data;
    delete[] weight;
  }

  //* 获取坐标[x, y]处的三通道值
  Vector3f getValue(const Vector2i &xy) const;

  //* 设置坐标[x, y]处的三通道值
  void addValue(const Vector2i &xy, const Vector3f &val, float weight);

  //* 以PNG格式保存该图片
  void savePNG(const char *filename) const;

  //* 以HDR格式保存该图片
  void saveHDR(const char *filename) const;

public:
  Vector2i size;
  static constexpr int channels = 3;

  friend std::shared_ptr<Image> loadImage(const char *filepath);
  friend std::string saveImage(std::string filepath,
                               std::shared_ptr<Image> image, bool overwrite);

private:
  float *data = nullptr;
  float *weight = nullptr;

  void normaliz();
};

//* 根据路径加载一张图片(PNG/JPG/HDR/EXR)
std::shared_ptr<Image> loadImage(const char *filepath);

//* 保存一张图片(PNG/JPG/HDR/EXR)
std::string saveImage(std::string filepath, std::shared_ptr<Image> image,
                      bool overwrite = false);
