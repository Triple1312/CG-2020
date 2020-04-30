//
// Created by Phili on 27/04/2020.
//

#ifndef ENGINE__PLATONICS_H_
#define ENGINE__PLATONICS_H_

#include "Figures3D.h"

class Platonic : public figures_3d::Figure {

};

class Cube : public Platonic {
 public:
  Cube();
};


 class Octahedron :  public Platonic {
  public:
  Octahedron();
};


 class Cone : public Platonic {
  public:
   Cone(const ini::Configuration &configuration, unsigned int i);
 };


 class Cylinder : public Platonic {
  public:
   Cylinder(const ini::Configuration &configuration, unsigned int i);
 };


 class Torus : public Platonic {
  public:
   Torus(const ini::Configuration &configuration, unsigned int i);
 };


 class Tetrahedron : public Platonic {
  public:
   Tetrahedron();
 };


 class Icosahedron : public Platonic {

   std::vector<std::vector<int>> ico_faces();

  public:

   Icosahedron();

 };


 class Dodecahedron : public Icosahedron {
  public:
   Dodecahedron(); // todo werkt niet

 };

 class Sphere : public Icosahedron {
  public:
   Sphere(int n);

 };

 //class MengerSponge() : public

#endif //ENGINE__PLATONICS_H_
