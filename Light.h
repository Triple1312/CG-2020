//
// Created by Phili on 27/04/2020.
//

#ifndef ENGINE__LIGHT_H_
#define ENGINE__LIGHT_H_

#include <list>

#include "easy_image.h"
#include "vector3d.h"
#include "ini_configuration.h"

class Light {
 public:

  Light();

  img::Color ambientLight;

  img::Color diffuseLight;

  img::Color specularLight;

  Vector3D ldVector;

  Vector3D location;

  double spotAngle;

};

class Lights3D : public std::list<Light> {
 public:
  Lights3D() = default;

  Lights3D(const ini::Configuration &config);

};

//class InfLight: public Light {
// public:
//  //de richting waarin het
//  //licht schijnt
//
//};
//
//class PointLight: public Light {
// public:
//  //de locatie van de puntbron
//
//   //de hoek van een spotlicht
//
//};


#endif //ENGINE__LIGHT_H_
