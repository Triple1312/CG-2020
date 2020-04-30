//
// Created by Phili on 27/04/2020.
//

#include "Platonics.h"

Cube::Cube() {
  points.push_back(Vector3D::point(1, -1, -1));
  points.push_back(Vector3D::point(-1, 1, -1));
  points.push_back(Vector3D::point(1, 1, 1));
  points.push_back(Vector3D::point(-1, -1, 1));
  points.push_back(Vector3D::point(1, 1, -1));
  points.push_back(Vector3D::point(-1, -1, -1));
  points.push_back(Vector3D::point(1, -1, 1));
  points.push_back(Vector3D::point(-1, 1, 1));
  faces.push_back(std::vector<int>{0, 5, 1, 4});
  faces.push_back(std::vector<int>{6, 3, 7, 2});
  faces.push_back(std::vector<int>{1, 5, 3, 7});
  faces.push_back(std::vector<int>{2, 6, 0, 4});
  faces.push_back(std::vector<int>{0, 5, 3, 6});
  faces.push_back(std::vector<int>{2, 7, 1, 4});
}


Octahedron::Octahedron() {
  points.push_back(Vector3D::point(1, 0, 0 ));
  points.push_back(Vector3D::point(0, 1, 0 ));
  points.push_back(Vector3D::point(-1, 0, 0 ));
  points.push_back(Vector3D::point(0, -1, 0 ));
  points.push_back(Vector3D::point(0, 0, -1 ));
  points.push_back(Vector3D::point(0, 0, 1 ));

  faces.push_back(std::vector<int>{1, 2, 6});
  faces.push_back(std::vector<int>{2, 3, 6});
  faces.push_back(std::vector<int>{3, 4, 6});
  faces.push_back(std::vector<int>{4, 1, 6});
  faces.push_back(std::vector<int>{2, 1, 5});
  faces.push_back(std::vector<int>{3, 2, 5});
  faces.push_back(std::vector<int>{4, 3, 5});
  faces.push_back(std::vector<int>{1, 4, 5});
}


Cone::Cone(const ini::Configuration &configuration, unsigned int i) {
  double height = configuration["Figure" + std::to_string(i)]["height"].as_double_or_die();
  int points_size = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
  auto top = Vector3D::point(0, 0, height);
  points.push_back(top);
  auto start_point = Vector3D::point(1, 0, 0); points.push_back(start_point);

  for (auto j = 1; j <= points_size; j++){
    auto point = Vector3D::point(std::cos(2 * M_PI * j  / points_size), std::sin(2 * M_PI * j  / points_size), 0);
    points.push_back(point);
    std::vector<int> face;
    face.push_back(j); face.push_back(j - 1); face.push_back(0);
    faces.push_back(face);
  }
  std::vector<int> last_face = {0, 1, int(points.size()-2)};
  faces.push_back(last_face);
}


Cylinder::Cylinder(const ini::Configuration &configuration, unsigned int i) {
  double height = configuration["Figure" + std::to_string(i)]["height"].as_double_or_die();
  int n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
  std::vector<int> firts_face = {0, 1};

  for (int i = 0; i < n; i++) {
    points.emplace_back(Vector3D::point(std::cos(2 * M_PI * i  / n),
                                            std::sin(2 * M_PI * i  / n), 0));
  }
  for (int i = 0; i < n; i++) {
    points.emplace_back(Vector3D::point(std::cos(2 * M_PI * i  / n),
                                            std::sin(2 * M_PI * i  / n), height));
  }
  for (int j = 0; j < n-1; j++) {
    faces.push_back(std::vector<int>{j, j+1, j+n+1, j+n});
  }
  faces.push_back(std::vector<int>{0, n-1, 2*n -1, n});
  std::vector<int> onder;
  std::vector<int> boven;
  for (int k = 0; k < n; k++) {
    onder.push_back(k);
    boven.push_back(n+k);
  }
  faces.push_back(onder);
  faces.push_back(boven);
}

Torus::Torus(const ini::Configuration &configuration, unsigned int i) {
  auto R = configuration["Figure" + std::to_string(i)]["R"].as_double_or_die();
  auto r = configuration["Figure" + std::to_string(i)]["r"].as_double_or_die();
  auto n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
  auto m = configuration["Figure" + std::to_string(i)]["m"].as_int_or_die();

  for (auto i = 0; i < n; i++){
    for (auto j = 0; j < m; j++){
      double u = 2 * i * M_PI / n;
      double v = 2 * j * M_PI / m;
      points.push_back(Vector3D::point(((R + r * std::cos(v)) * std::cos(u)),
                                           ((R + r * std::cos(v)) * std::sin(u)),
                                           r * std::sin(v)));
    }
  }
  for (auto i = 0; i < n; i++){
    for (auto j = 0; j < m; j++){
      faces.push_back(std::vector<int>{i*m+j, ((i+1)%n)*m+j,((i+1)%n)*m+((j+1)%m), i*m+((j+1)%m)});
    }
  }
}

Tetrahedron::Tetrahedron() {
  points.push_back(Vector3D::point(1, -1, -1));
  points.push_back(Vector3D::point(-1, 1, -1));
  points.push_back(Vector3D::point(1, 1, 1));
  points.push_back(Vector3D::point(-1, -1, 1));
  faces.push_back(std::vector<int>{0, 1, 2});
  faces.push_back(std::vector<int>{0, 2, 3});
  faces.push_back(std::vector<int>{0, 1, 3});
  faces.push_back(std::vector<int>{1, 2, 3});
}


inline std::vector<std::vector<int>> Icosahedron::ico_faces() {
  std::vector<std::vector<int>> faces;
  faces.push_back(std::vector<int>{0, 1, 2}); faces.push_back(std::vector<int>{0, 2, 3});
  faces.push_back(std::vector<int>{0, 3, 4}); faces.push_back(std::vector<int>{0, 4, 5});
  faces.push_back(std::vector<int>{0, 5, 1}); faces.push_back(std::vector<int>{1, 6, 2});
  faces.push_back(std::vector<int>{1, 10, 6}); faces.push_back(std::vector<int>{2, 6, 7});
  faces.push_back(std::vector<int>{2, 7, 3}); faces.push_back(std::vector<int>{3, 7, 8});
  faces.push_back(std::vector<int>{3, 8, 4}); faces.push_back(std::vector<int>{4, 8, 9});
  faces.push_back(std::vector<int>{4, 9, 5}); faces.push_back(std::vector<int>{5, 9, 10});
  faces.push_back(std::vector<int>{5, 10, 1}); faces.push_back(std::vector<int>{11, 6, 10});
  faces.push_back(std::vector<int>{11, 7, 6}); faces.push_back(std::vector<int>{11, 8, 7});
  faces.push_back(std::vector<int>{11, 9, 8}); faces.push_back(std::vector<int>{11, 10, 9});
  return faces;
}

Icosahedron::Icosahedron() {
  points.push_back(Vector3D::point(0, 0, sqrt(5) / 2));
  for (auto i = 2; i <= 6; i++) {
    points.push_back(Vector3D::point(std::cos((i-2) * 2 * M_PI / 5),
                                         std::sin((i-2) * 2 * M_PI / 5),
                                         0.5));
  }
  for (int j = 7; j <= 11 ; j++) {
    points.push_back(Vector3D::point(std::cos( M_PI / 5 + (j-7) * 2 * M_PI / 5),
                                         std::sin( M_PI / 5 + (j-7) * 2 * M_PI / 5),
                                         -0.5));
  }
  points.push_back(Vector3D::point(0, 0, - sqrt(5) / 2));

  faces = ico_faces();
}


Dodecahedron::Dodecahedron() {
  figures_3d::Figure ico = Icosahedron();

  for (auto i : ico.faces){
    
    this->points.push_back(Vector3D::point(( ico.points[i[0]].x + ico.points[i[1]].x + ico.points[i[2]].x) / 3,
                                         ( ico.points[i[0]].y + ico.points[i[1]].y + ico.points[i[2]].y) / 3,
                                         ( ico.points[i[0]].z + ico.points[i[1]].z + ico.points[i[2]].z) / 3) );
  }

  for (int j = 0; j < this->points.size(); j++) {
    for (int k = j+1; k < this->points.size(); k++) {
      if ( sqrt(pow(this->points[j].x + this->points[k].x, 2) +
          pow(this->points[j].y + this->points[k].y, 2) +
          pow(this->points[j].z + this->points[k].z, 2)) > 1.5){
        this->faces.push_back(std::vector<int>{j, k}); //
      }
    }
  }

}




Sphere::Sphere(int n) {

  for (auto j = 0; j < n; j++) {
    auto faces_b_s = this->faces.size(); // begin size of faces
    for (int i = 0; i < faces_b_s; i++) {
      this->points.push_back((this->points[this->faces[0][0]] + this->points[this->faces[0][1]]) / 2);
      this->points.push_back((this->points[this->faces[0][0]] + this->points[this->faces[0][2]]) / 2);
      this->points.push_back((this->points[this->faces[0][1]] + this->points[this->faces[0][2]]) / 2);

      this->faces.push_back(
          std::vector<int>{this->faces[0][0], int(this->points.size()) - 3, int(this->points.size()) - 2});
      this->faces.push_back(
          std::vector<int>{this->faces[0][1], int(this->points.size()) - 3, int(this->points.size()) - 1});
      this->faces.push_back(
          std::vector<int>{this->faces[0][2], int(this->points.size()) - 2, int(this->points.size()) - 1});
      this->faces.push_back(std::vector<int>{int(this->points.size()) - 1, int(this->points.size()) - 3,
                                           int(this->points.size()) - 2});
      this->faces.erase(this->faces.begin());
    }
  }

  for (int k = 0; k < this->points.size(); k++) {
    double r = sqrt(pow(this->points[k].x, 2) +
        pow(this->points[k].y, 2) +
        pow(this->points[k].z, 2));
    this->points[k] = this->points[k] / r;
  }
  
}
