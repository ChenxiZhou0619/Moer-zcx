#include <ResourceLayer/Image.h>
#include <regex>
//* Apply postprocesses on output image

//* Linear
std::shared_ptr<Image> filter(std::shared_ptr<Image> image) {
  std::shared_ptr<Image> filtered = std::make_shared<Image>(image->size);

  auto filter_sample = [&](int x, int y, int kernel_size) {
    int len = kernel_size / 2;
    Vector3f value(.0f);
    for (int dx = -len; dx <= len; ++dx) {
      for (int dy = -len; dy <= len; ++dy) {
        value += image->getValue({x + dx, y + dy});
      }
    }
    float weight = (1 + 2 * len) * (1 + 2 * len);
    return value / weight;
  };

  auto filter_sample_gaussian = [&](int x, int y, int kernel_size) {
    int len = kernel_size / 2;
    Vector3f value(.0f);
    float weight = .0f;
    float determinant = 1.f / (2.f * PI * len);
    for (int dx = -len; dx <= len; ++dx) {
      for (int dy = -len; dy <= len; ++dy) {
        float dist2 = dx * dx + dy * dy;
        float w = determinant * std::exp(-.5f * dist2 / (len * len));
        value += image->getValue({x + dx, y + dy}) * w;
        weight += w;
      }
    }
    return value / weight;
  };

  auto filter_sample_bilateral = [&](int x, int y, int kernel_size,
                                     float sigma_r) {
    int len = kernel_size / 2;
    Vector3f value(.0f);
    Vector3f weight(.0f);

    float determinant_s = 1.f / (2.f * PI * len);
    float determinant_r = 1.f / (2.f * PI * sigma_r);

    Vector3f rgb = image->getValue({x, y});

    for (int dx = -len; dx <= len; ++dx) {
      for (int dy = -len; dy <= len; ++dy) {
        float dist2 = dx * dx + dy * dy;
        float w_s = determinant_s * std::exp(-.5f * dist2 / (len * len));
        Vector3f rgb_new = image->getValue({x + dx, y + dy});
        Vector3f w_r =
            determinant_r *
            Vector3f{std::exp(-.5f * std::pow((rgb[0] - rgb_new[0]), 2.f)),
                     std::exp(-.5f * std::pow((rgb[1] - rgb_new[1]), 2.f)),
                     std::exp(-.5f * std::pow((rgb[2] - rgb_new[2]), 2.f))};
        w_r /= (sigma_r * sigma_r);

        value += rgb_new * w_s * w_r;
        weight += w_s * w_r;
      }
    }
    return value / weight;
  };

  //* Filter each pixel
  for (int x = 0; x < image->size[0]; ++x) {
    for (int y = 0; y < image->size[1]; ++y) {
      Vector3f v = filter_sample_gaussian(x, y, 3);
      // Vector3f v = filter_sample_bilateral(x, y, 5, 0.2);
      filtered->addValue({x, y}, v, 1.f);
    }
  }
  return filtered;
}

int main(int argc, char **argv) {
  if (argv[1] == nullptr) {
    printf("Image should be provided!\n");
    exit(1);
  }

  //*----- Load Image -----
  auto image = loadImage(argv[1]);

  //*----- Output filename -----
  std::string path(argv[1]);
  std::regex rgx("(.*)(\\.)(.*)");
  auto dest = "$1-" + std::string("filter") + ".$3";
  auto res = std::regex_replace(path, rgx, dest);

  //*----- Filter -----
  auto image_after_filter = filter(image);

  printf("%s\n", res.c_str());

  //*----- Save Output -----
  saveImage(res, image_after_filter, true);
}