#include "CheckboardTexture.h"
#include <ResourceLayer/Factory.h>
CheckboardTexture::CheckboardTexture(const Json &json) {
  size = fetchOptional(json, "size", 0.01f);
  data1 = fetchOptional(json, "data1", Spectrum(.1f));
  data2 = fetchOptional(json, "data2", Spectrum(.7f));
}

REGISTER_CLASS(CheckboardTexture, "checkboardTex")