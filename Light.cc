//
// Created by Phili on 27/04/2020.
//

#include "Light.h"


Lights3D::Lights3D(const ini::Configuration &config) {
  for (auto i = 0; i < config["General"]["nrLights"].as_int_or_default(0); i++) {
    ini::DoubleTuple am = config["Light"+ std::to_string(i)]["ambientLight"].as_double_tuple_or_default(ini::DoubleTuple(3, 0));
    ini::DoubleTuple dif = config["Light"+ std::to_string(i)]["diffuseLight"].as_double_tuple_or_default(ini::DoubleTuple(3, 0));
    ini::DoubleTuple dir = config["Light"+ std::to_string(i)]["direction"].as_double_tuple_or_default(ini::DoubleTuple(3, 0));
    Light light; light.ambientLight = img::Color(am[0], am[1], am[2]);
    light.diffuseLight = img::Color(dif[0], dif[1], dif[2]);
    light.ldVector = Vector3D::vector(dir[0], dir[1], dir[2]);
    this->emplace_back(light);
  }
}

Light::Light() {
  this->ambientLight = img::Color(255, 255, 255);
}
