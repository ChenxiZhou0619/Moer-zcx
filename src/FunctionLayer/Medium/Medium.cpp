#include "Medium.h"
#include "./Phase/HGPhase.h"
#include "./Phase/IsotropicPhase.h"
Medium::Medium(const Json &json) {
  if (json.contains("g")) {
    float g = fetchRequired<float>(json, "g");
    phase = std::make_shared<HenyeyGrennsteinPhase>(g);
  } else {
    phase = std::make_shared<IsotropicPhase>();
  }
}