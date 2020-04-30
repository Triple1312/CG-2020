//
// Created by Phili on 3/4/2020.
//

#include "Figures3D.h"

#include "Platonics.h"

namespace {
    inline void rot_U_axis(const double angle, std::tuple<double, double, double> &H_axis, std::tuple<double, double, double> &L_axis) {
      std::tuple<double, double, double> H_new(std::get<0>(H_axis) * cos(angle) + std::get<0>(L_axis) * sin(angle),
                                               std::get<1>(H_axis) * cos(angle) + std::get<1>(L_axis) * sin(angle),
                                               std::get<2>(H_axis) * cos(angle) + std::get<2>(L_axis) * sin(angle));
      std::tuple<double, double, double> L_new(std::get<0>(L_axis) * cos(angle) - std::get<0>(H_axis) * sin(angle),
                                               std::get<1>(L_axis) * cos(angle) - std::get<1>(H_axis) * sin(angle),
                                               std::get<2>(L_axis) * cos(angle) - std::get<2>(H_axis) * sin(angle));
      H_axis = H_new; L_axis = L_new;
    }

    inline void rot_L_axis(const double angle, std::tuple<double, double, double> &H_axis, std::tuple<double, double, double> &U_axis) {
      std::tuple<double, double, double> H_new(std::get<0>(H_axis) * cos(angle) + std::get<0>(U_axis) * sin(angle),
                                               std::get<1>(H_axis) * cos(angle) + std::get<1>(U_axis) * sin(angle),
                                               std::get<2>(H_axis) * cos(angle) + std::get<2>(U_axis) * sin(angle));
      std::tuple<double, double, double> U_new(std::get<0>(U_axis) * cos(angle) - std::get<0>(H_axis) * sin(angle),
                                               std::get<1>(U_axis) * cos(angle) - std::get<1>(H_axis) * sin(angle),
                                               std::get<2>(U_axis) * cos(angle) - std::get<2>(H_axis) * sin(angle));
      H_axis = H_new; U_axis = U_new;
    }

    inline void rot_H_axis(const double angle, std::tuple<double, double, double> &L_axis, std::tuple<double, double, double> &U_axis) {
      std::tuple<double, double, double> L_new(std::get<0>(L_axis) * cos(angle) - std::get<0>(U_axis) * sin(angle),
                                               std::get<1>(L_axis) * cos(angle) - std::get<1>(U_axis) * sin(angle),
                                               std::get<2>(L_axis) * cos(angle) - std::get<2>(U_axis) * sin(angle));
      std::tuple<double, double, double> U_new(std::get<0>(U_axis) * cos(angle) + std::get<0>(L_axis) * sin(angle),
                                               std::get<1>(U_axis) * cos(angle) + std::get<1>(L_axis) * sin(angle),
                                               std::get<2>(U_axis) * cos(angle) + std::get<2>(L_axis) * sin(angle));
      L_axis = L_new; U_axis = U_new;
    }

    inline figures_3d::Figure line_drawing(const ini::Configuration &configuration, unsigned int i){
        figures_3d::Figure figure;
        for (auto j=0; j < configuration["Figure" + std::to_string(i)]["nrPoints"].as_int_or_die() ; j++){
            ini::DoubleTuple point = configuration["Figure" + std::to_string(i)]["point" + std::to_string(j)].as_double_tuple_or_die();
            //Vector3D newPoint = Vector3D::point(point[0], point[1], point[2]);
            figure.points.push_back(Vector3D::point(point[0], point[1], point[2]));
        }

        for (auto k = 0; k < configuration["Figure" + std::to_string(i)]["nrLines"].as_int_or_die(); k++ ) {
            ini::IntTuple line = configuration["Figure" + std::to_string(i)]["line" + std::to_string(k)].as_int_tuple_or_die();
            //Vector3D newPoint = Vector3D::point(point[0], point[1], point[2]);
            std::vector<int> line_n; line_n.push_back(line[0]); line_n.push_back(line[1]);
            figure.faces.push_back(line_n);
        }
        return figure;
    }

    inline Matrix getMatrix_rotate_x(const double &angle){
        Matrix rotate_matrix;
        rotate_matrix( 2, 2) = cos(angle);
        rotate_matrix( 2, 3) = sin(angle);
        rotate_matrix( 3, 2) = -sin(angle);
        rotate_matrix( 3, 3) = cos(angle);
        return rotate_matrix;
    }

    inline Matrix getMatrix_rotate_y(const double &angle){
        Matrix rotate_matrix;
        rotate_matrix( 1, 1) = cos(angle);
        rotate_matrix( 3, 1) = sin(angle);
        rotate_matrix( 1, 3) = -sin(angle);
        rotate_matrix( 3, 3) = cos(angle);
        return rotate_matrix;
    }

    inline Matrix getMatrix_rotate_z(const double &angle){
        Matrix rotate_matrix;
        rotate_matrix( 1, 1) = cos(angle);
        rotate_matrix( 1, 2) = sin(angle);
        rotate_matrix( 2, 1) = -sin(angle);
        rotate_matrix( 2, 2) = cos(angle);
        return rotate_matrix;
    }

    inline Matrix getMatrix_scale(const double& scale ) {
      Matrix scaler;
      scaler( 1, 1 ) = scale;
      scaler( 2, 2 ) = scale;
      scaler( 3, 3 ) = scale;
      scaler( 4, 4 ) = 1;
      return scaler;
    }

    inline Matrix getMatrix_translate(const double &x, const double &y, const double &z){
        Matrix translate_matrix;
        translate_matrix( 4, 1) = x;
        translate_matrix( 4, 2) = y;
        translate_matrix( 4, 3) = z;
        return translate_matrix;
    }

    inline Matrix getMatrix_translate(const Vector3D& vector){
      Matrix translate_matrix;
      translate_matrix( 4, 1) = vector.x;
      translate_matrix( 4, 2) = vector.y;
      translate_matrix( 4, 3) = vector.z;
      return translate_matrix;
    }

    inline Matrix eye_matrix( const double &theta, const double &phi, const double &r){
        Matrix eye_matrix;
        eye_matrix(1, 1) = -sin(theta);
        eye_matrix(1, 2) = -cos(theta) * cos(phi);
        eye_matrix(1, 3) = cos(theta) * sin(phi);
        eye_matrix(2,1) = cos(theta);
        eye_matrix(2, 2) = -sin(theta) * cos(phi);
        eye_matrix(2, 3) = sin(theta) * sin(phi);
        eye_matrix(3, 2) = sin(phi);
        eye_matrix(3, 3) = cos(phi);
        eye_matrix(4, 3) = -r;
        return eye_matrix;
    }

    inline Matrix eye_clip_matrix(const double &theta, const double &phi, Vector3D &eye){
      return Matrix(getMatrix_translate(-eye) * getMatrix_rotate_z(-(theta + M_PI / 2)) * getMatrix_rotate_x(- phi));
    }

    ///engine crasht als nrIt == 0
figures_3d::Figures3D menger_sponge(int nrIt, figures_3d::Figure fig = Cube() ) {
  figures_3d::Figures3D figs;
  std::vector<double> ops;
  ops.emplace_back( 2.0/3);
  ops.emplace_back(-2.0/3);
  ops.emplace_back(0);
  nrIt -= 1;

  for (auto &i : ops) {
    for (auto &j : ops) {
      for (auto &k : ops) {
        if (i*j != 0|| j*k != 0 || i*k != 0){
          auto cube_n = fig;
          cube_n.scale(1.0/3);
          cube_n.translate(i, j, k);
          //cube_n.ambientReflection = color; // todo check ook nu werkt
          figs.emplace_back(cube_n);
        }
      }
    }
  }
  if (nrIt != 0) {
    figures_3d::Figures3D temps;
    for (auto &f: figs){
      auto temps2 = menger_sponge(nrIt -1, f);
      temps.insert(temps.end(), temps2.begin(), temps2.end());
    }
    figs = temps;
  }
  return figs;
}

} // namespace

figures_3d::Figures3D figures_3d::fractal(int nrIt, const double& scale, Figure fig = Cube()) {

  figures_3d::Figures3D temps;
  Matrix scaler = getMatrix_scale(1 / scale);
  if (nrIt == 0) {return {fig};}
    for (auto j : fig.points) {
      Vector3D Pai = j * scaler;
      Matrix Mt = getMatrix_translate(j - Pai);
      if (1 < nrIt) {
        figures_3d::Figures3D new_figs = fractal(nrIt-1, scale, fig);
        for (auto &m : new_figs) {
          for (auto &l : m.points) {
            l *= scaler;
            l *= Mt;
          }
          temps.push_back(m);
        }
      }
      else {
        auto fig_n = fig;
        for ( auto &l : fig_n.points) {
          l *= scaler;
          l *= Mt;
        }
        temps.push_back(fig_n);
      }
    }
    return temps;
}

figures_3d::Figure figures_3d::calc_fig(LParser::LSystem3D & l_sys) {
  figures_3d::Figure fig;

  fig.points.emplace_back(Vector3D::point(0.0, 0.0, 0.0));
  std::string full_function = L_system::getString(l_sys);

  std::tuple<double, double, double> H(1, 0, 0);
  std::tuple<double, double, double> L(0, 1, 0);
  std::tuple<double, double, double> U(0, 0, 1);
  std::stack<std::tuple<double, double, double, std::tuple<double, double, double>, std::tuple<double, double, double>, std::tuple<double, double, double>>> stack;

  double angle = l_sys.get_angle() * M_PI /180;
  double curr_x(0.0), curr_y(0.0), curr_z(0.0);
  int test = 0;

  int counter = 0;
  for (char& i : full_function) {

    if (i == '+') {
      rot_U_axis(angle, H, L);
    } else if (i == '-') {
      rot_U_axis(-angle, H, L);
    } else if (i == '^') {
      rot_L_axis(angle, H, U);
    } else if (i == '&') {
      rot_L_axis(-angle, H, U);
    } else if (i == '\\') {
      rot_H_axis(angle, L, U);
    } else if (i == '/') {
      rot_H_axis(-angle, L, U);
    } else if (i == '(') {
      stack.push(std::tuple<double,
                            double,
                            double,
                            std::tuple<double, double, double>,
                            std::tuple<double, double, double>,
                            std::tuple<double, double, double>>
                     (curr_x, curr_y, curr_z, H, L, U));

    } else if (i == ')') {
      std::tie(curr_x, curr_y, curr_z, H, L, U) = stack.top();
      fig.points.emplace_back(Vector3D::point(curr_x, curr_y, curr_z));
      stack.pop();
    } else {
      test += 1;
      std::cout << test << i << std::endl;
      if (l_sys.draw(i)) {
        counter += 1;
        curr_x += std::get<0>(H);
        curr_y += std::get<1>(H);
        curr_z += std::get<2>(H);
        fig.points.emplace_back(Vector3D::point(curr_x, curr_y, curr_z));
        fig.faces.emplace_back(figures_3d::Face{(int) fig.points.size() - 2, (int) fig.points.size() - 1});
      }
    }
    std::cout << test << std::endl;
    std::cout << counter << std::endl;
  }
  return fig;
}

void figures_3d::Figure::scale(const double &scale) {
    Matrix scaler;
    scaler( 1, 1 ) = scale;
    scaler( 2, 2 ) = scale;
    scaler( 3, 3 ) = scale;
    scaler( 4, 4 ) = 1;

    for ( auto &point : this->points) {
        point *= scaler;
    }
}


void figures_3d::Figure::rotateX(const double &angle) {
    //Matrix rot_matrix = getMatrix_rotate_x(angle);
    for ( auto &point : this->points) {
        point *= getMatrix_rotate_x(angle);
    }
}

void figures_3d::Figure::rotateY(const double &angle) {
    auto rot_matrix = getMatrix_rotate_y(angle);
    for ( auto &point : this->points) {
        point *= rot_matrix;
    }
}

void figures_3d::Figure::rotateZ(const double &angle) {
    auto rot_matrix = getMatrix_rotate_z(angle);
    for ( auto &point : this->points) {
        point *= rot_matrix;
    }
}

void figures_3d::Figure::translate(const double &x, const double &y, const double &z) {
    auto trans_matrix = getMatrix_translate( x, y, z );
    for ( auto &point : this->points) {
        point *= trans_matrix;
    }
}

lines_2d::Lines2D figures_3d::Figure::project_figure(double d) {
    lines_2d::Lines2D lines;
    for (auto k = 0; k < this->faces.size(); k++){
        for (auto i = 0; i < this->faces[k].size()-1; i ++){
            Vector3D &point1 = this->points[this->faces[k][i]];
            Vector3D &point2 = this->points[this->faces[k][i+1]];
            lines.push_back(lines_2d::Line2D(
                    lines_2d::Point2D(- point1.x * d / point1.z , - point1.y * d / point1.z),
                    lines_2d::Point2D(- point2.x * d / point2.z , - point2.y * d / point2.z),
                    this->ambientReflection, point1.z, point2.z ));

        }
        Vector3D &point1 = this->points[this->faces[k][0]];
        Vector3D &point2 = this->points[this->faces[k][this->faces[k].size()-1]];
        lines.push_back(lines_2d::Line2D(
            lines_2d::Point2D(- point1.x * d / point1.z , - point1.y * d / point1.z),
            lines_2d::Point2D(- point2.x * d / point2.z , - point2.y * d / point2.z),
            this->ambientReflection, point1.z, point2.z ));
    }
    return lines;
}

void figures_3d::Figure::to_eye( const double &theta, const double &phi, const double &r) {
    for (auto &i : this->points){
        i *= eye_matrix(theta, phi, r);
    }
}

void figures_3d::Figure::to_eye_clip(const double &theta, const double &phi, Vector3D &eye) {
  for (auto &i :this->points) {
    i *= eye_clip_matrix(theta, phi, eye );
  }
}

img::Color figures_3d::Figure::calcColor(Lights3D &lights) {
  img::Color color;
  for(auto &i : lights){
    color.red = color.red + ambientReflection.red * i.ambientLight.red;
    color.blue = color.blue + ambientReflection.blue * i.ambientLight.blue;
    color.green = color.green + ambientReflection.green * i.ambientLight.green;
  }

  return color;
}



figures_3d::Figures3D figures_3d::Wireframe(const ini::Configuration &configuration) {
    figures_3d::Figures3D figs;
    int blub = configuration["General"]["nrFigures"].as_int_or_die();

    ini::DoubleTuple eye = configuration["General"]["eye"].as_double_tuple_or_die();

    for (auto i = 0; i < configuration["General"]["nrFigures"].as_int_or_die(); i++){
      if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalCube" ||
          configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalDodecahedron" ||
          configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalTetrahedron" ||
          configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalIcosahedron" ||
          configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "MengerSponge"){
        int nrIt = configuration["Figure" + std::to_string(i)]["nrIterations"].as_int_or_die();
        figures_3d::Figures3D temp_figs;
        if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalCube") {
          temp_figs = fractal(nrIt,configuration["Figure"+ std::to_string(i)]["fractalScale"].as_double_or_die(), Cube());
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() =="FractalTetrahedron"){
          temp_figs = fractal(nrIt,configuration["Figure"+ std::to_string(i)]["fractalScale"].as_double_or_die(), Tetrahedron());
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalIcosahedron"){
          temp_figs = fractal(nrIt,configuration["Figure"+ std::to_string(i)]["fractalScale"].as_double_or_die(), Icosahedron());
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalDodecahedron"){
          temp_figs = fractal(nrIt,configuration["Figure"+ std::to_string(i)]["fractalScale"].as_double_or_die(), Dodecahedron());
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "MengerSponge"){
          auto color = configuration["Figure" + std::to_string(i)]["color"].as_double_tuple_or_die();
          auto Color = img::Color(color[0]* 255, color[1] * 255, color[2] *255);
          auto sponge_fig = Cube(); sponge_fig.ambientReflection = Color;
          temp_figs = menger_sponge(nrIt ,sponge_fig);
        }
        for (auto k : temp_figs) {
          ini::DoubleTuple am_refl = configuration["Figure" + std::to_string(i)]["ambientReflection"].as_double_tuple_or_default(ini::DoubleTuple(3, 0));
          k.ambientReflection = img::Color(am_refl[0] * 255, am_refl[1] *255, am_refl[2]*255);
          ini::DoubleTuple fig_color = configuration["Figure" + std::to_string(i)]["color"].as_double_tuple_or_default(ini::DoubleTuple(3, 0));
          k.ambientReflection = img::Color(k.ambientReflection.red +fig_color[0] * 255, k.ambientReflection.green + fig_color[1] * 255, k.ambientReflection.blue + fig_color[2] * 255);
          double r = std::pow(std::pow(eye[0], 2.0) + std::pow(eye[1], 2.0) + std::pow(eye[2], 2.0), 0.5);
          ini::DoubleTuple center = configuration["Figure" + std::to_string(i)]["center"].as_double_tuple_or_die();

          //k.scale(configuration["Figure" + std::to_string(i)]["scale"].as_double_or_die());

          k.rotateX(M_PI / 180 * configuration["Figure" + std::to_string(i)]["rotateX"].as_double_or_die());
          k.rotateY(M_PI / 180 * configuration["Figure" + std::to_string(i)]["rotateY"].as_double_or_die());
          k.rotateZ(M_PI / 180 * configuration["Figure" + std::to_string(i)]["rotateZ"].as_double_or_die());
          k.translate(center[0], center[1], center[2]);
          if (!configuration["General"]["clipping"].as_bool_or_default(false)) {
            k.to_eye(std::atan2(eye[1], eye[0]), std::acos(eye[2] / r), r);
          }
          //k.to_eye(std::atan2(eye[1], eye[0]), std::acos(eye[2] / r), r);
          figs.push_back(k);
        }
      }
      else {

        figures_3d::Figure figure;  // todo kan veel beter door direct to initializen in de vector
        if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "LineDrawing") {
          figure = line_drawing(configuration, i);  // i is the number of the figure
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Cube") {
          figure = Cube();
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Cone") {
          figure = Cone(configuration, i);
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Cylinder") {
          figure = Cylinder(configuration, i);
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Torus") {
          figure = Torus(configuration, i);
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Tetrahedron") {
          figure = Tetrahedron();
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Icosahedron") {
          figure = Icosahedron();
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Octahedron") {
          figure = Octahedron();
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Dodecahedron") {
          figure = Dodecahedron();
        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Sphere") {
          figure = Sphere(configuration["Figure" + std::to_string(i)]["n"].as_int_or_die());
        }
//        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "FractalCube") {
//          figs = fractal_cube(configuration["Figure" + std::to_string(i)]["nrIterations"].as_int_or_die(),
//                              configuration["Figure"
//                                  + std::to_string(i)]["fractalScale"].as_int_or_die()); //todo color meegeven
//          break;
//        }
        else if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "3DLSystem") {
          std::ifstream file(configuration["Figure" + std::to_string(i)]["inputfile"].as_string_or_die());
          LParser::LSystem3D l_sys(file);
          ini::DoubleTuple lineColor = configuration["Figure" + std::to_string(i)]["color"].as_double_tuple_or_default(ini::DoubleTuple(3, 0));
          img::Color line_c(lineColor[0]*255, lineColor[1]*255, lineColor[2]*255);
          figure = figures_3d::calc_fig(l_sys);
          //figure.color = line_c;

        }

        ini::DoubleTuple fig_color = configuration["Figure" + std::to_string(i)]["color"].as_double_tuple_or_default(ini::DoubleTuple(3, 0));
        ini::DoubleTuple am_refl = configuration["Figure" + std::to_string(i)]["ambientReflection"].as_double_tuple_or_default(ini::DoubleTuple(3, 0));
        figure.ambientReflection = img::Color(am_refl[0]*255, am_refl[1]*255, am_refl[2]*255);
        figure.ambientReflection = img::Color(figure.ambientReflection.red + fig_color[0] * 255, figure.ambientReflection.green +fig_color[1] * 255, figure.ambientReflection.blue +fig_color[2] * 255);
        double r = std::pow(std::pow(eye[0], 2.0) + std::pow(eye[1], 2.0) + std::pow(eye[2], 2.0), 0.5);
        ini::DoubleTuple center = configuration["Figure" + std::to_string(i)]["center"].as_double_tuple_or_die();

        figure.scale(configuration["Figure" + std::to_string(i)]["scale"].as_double_or_die());

        figure.rotateX(M_PI / 180 * configuration["Figure" + std::to_string(i)]["rotateX"].as_double_or_die());
        figure.rotateY(M_PI / 180 * configuration["Figure" + std::to_string(i)]["rotateY"].as_double_or_die());
        figure.rotateZ(M_PI / 180 * configuration["Figure" + std::to_string(i)]["rotateZ"].as_double_or_die());
        figure.translate(center[0], center[1], center[2]);
        if (!configuration["General"]["clipping"].as_bool_or_default(false)) {
          figure.to_eye(std::atan2(eye[1], eye[0]), std::acos(eye[2] / r), r);
        }

        figs.push_back(figure);
      }
    }
    return figs;
}

std::vector<figures_3d::Face> figures_3d::triangulate(figures_3d::Face & face) {
  std::vector<Face> new_faces;
  for( auto i = 1; i < face.size() - 1; i++ ){
    Face new_face; new_face.push_back(face[0]);
    new_face.push_back(face[i]); new_face.push_back(face[i+1]);
    new_faces.push_back(new_face);
  }
  return new_faces;
}

std::vector<figures_3d::Face> figures_3d::triangulate(std::vector<figures_3d::Face> &faces ) {
  std::vector<Face> new_faces;
  for(auto i : faces){
      for (auto j : triangulate(i)){
        new_faces.emplace_back(j);
      }
  }
  return new_faces;
}

void figures_3d::draw_triangle(const Vector3D &A,
                               const Vector3D &B,
                               const Vector3D &C,
                               double dx,
                               double dy,
                               img::EasyImage &image,
                               zbuffer::ZBuffer &buffer,
                               img::Color& color,
                               double d) {

  lines_2d::Point2D pa(d * A.x / (- A.z) + dx, d * A.y / (- A.z) + dy);
  lines_2d::Point2D pb(d * B.x / (- B.z) + dx, d * B.y / (- B.z) + dy);
  lines_2d::Point2D pc(d * C.x / (- C.z) + dx, d * C.y / (- C.z) + dy);
  std::vector<double> y; y.push_back(pa.y);
  y.push_back(pb.y); y.push_back(pc.y);

  for (int yi = tools::d2i(min(y)+0.5); yi <= tools::d2i(max(y) -0.5); yi++) {
    std::vector<double> xl = {std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity(),
                   std::numeric_limits<double>::infinity()};
    std::vector<double> xr = { - std::numeric_limits<double>::infinity(),
                    - std::numeric_limits<double>::infinity(),
                    - std::numeric_limits<double>::infinity()};
    if ( (yi - pb.y)*(yi - pa.y) <= 0) {
      xl[0] = xr[0] = pa.x + (pb.x - pa.x) * (yi - pa.y) / (pb.y - pa.y);
    }
    if ( (yi - pc.y)*(yi - pa.y) <= 0) {
      xl[1] = xr[1] = pa.x + (pc.x - pa.x) * (yi - pa.y) / (pc.y - pa.y);
    }
    if ( (yi - pb.y)*(yi - pc.y) <= 0) {
      xl[2] = xr[2] = pc.x + (pb.x - pc.x) * (yi - pc.y) / (pb.y - pc.y);
    }
    int xl_ = tools::d2i(min(xl) +0.5);
    int xr_ = tools::d2i(max(xr) -0.5);
    Vector3D u = B - A;
    Vector3D v = C - A;
    Vector3D w = Vector3D::cross(u, v);
    double k = w.x * A.x + w.y * A.y + w.z * A.z;
    double dzdx = - w.x /d /k;
    double dzdy = - w.y /d /k;
    double xg = (pa.x + pb.x + pc.x) /3;
    double yg = (pa.y + pb.y + pc.y) /3;
    double ozg = 1 / (3 * A.z) + 1/ (3 * B.z) + 1/ (3 * C.z); // 1/Zg
    for (int xi = xl_; xi <= xr_; xi++) {
      double oz = 1.0001 * ozg + (xi - xg) * dzdx + (yi - yg) * dzdy;
      if (buffer(xi, yi) > oz){
        image(xi, yi) = color;
        buffer(xi, yi) = oz;
      }
    }
  }
}

void figures_3d::draw_triangle_light(const Vector3D &A,
                                     const Vector3D &B,
                                     const Vector3D &C,
                                     double dx,
                                     double dy,
                                     double d,
                                     img::EasyImage &image,
                                     zbuffer::ZBuffer &buffer,
                                     img::Color ambientReflection,
                                     img::Color diffuseReflection,
                                     img::Color specularReflection, double reflectionCoeff,
                                     Lights3D& lights) {

  lines_2d::Point2D pa(d * A.x / (- A.z) + dx, d * A.y / (- A.z) + dy);
  lines_2d::Point2D pb(d * B.x / (- B.z) + dx, d * B.y / (- B.z) + dy);
  lines_2d::Point2D pc(d * C.x / (- C.z) + dx, d * C.y / (- C.z) + dy);
  std::vector<double> y; y.push_back(pa.y);
  y.push_back(pb.y); y.push_back(pc.y);

  for (int yi = tools::d2i(min(y)+0.5); yi <= tools::d2i(max(y) -0.5); yi++) {
    std::vector<double> xl = {std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::infinity(),
                              std::numeric_limits<double>::infinity()};
    std::vector<double> xr = { - std::numeric_limits<double>::infinity(),
                               - std::numeric_limits<double>::infinity(),
                               - std::numeric_limits<double>::infinity()};
    if ( (yi - pb.y)*(yi - pa.y) <= 0) {
      xl[0] = xr[0] = pa.x + (pb.x - pa.x) * (yi - pa.y) / (pb.y - pa.y);
    }
    if ( (yi - pc.y)*(yi - pa.y) <= 0) {
      xl[1] = xr[1] = pa.x + (pc.x - pa.x) * (yi - pa.y) / (pc.y - pa.y);
    }
    if ( (yi - pb.y)*(yi - pc.y) <= 0) {
      xl[2] = xr[2] = pc.x + (pb.x - pc.x) * (yi - pc.y) / (pb.y - pc.y);
    }
    int xl_ = tools::d2i(min(xl) +0.5); // x left
    int xr_ = tools::d2i(max(xr) -0.5); // x right
    Vector3D u = B - A;
    Vector3D v = C - A;
    Vector3D w = Vector3D::cross(u, v);
    Vector3D n = Vector3D::normalise(w);

    double k = w.x * A.x + w.y * A.y + w.z * A.z;
    double dzdx = - w.x /d /k;
    double dzdy = - w.y /d /k;
    double xg = (pa.x + pb.x + pc.x) /3;
    double yg = (pa.y + pb.y + pc.y) /3;
    double ozg = 1 / (3 * A.z) + 1/ (3 * B.z) + 1/ (3 * C.z); // 1/Zg
    for (int xi = xl_; xi <= xr_; xi++) {
      double oz = 1.0001 * ozg + (xi - xg) * dzdx + (yi - yg) * dzdy;
      if (buffer(xi, yi) > oz){
        img::Color color;
        for ( auto &i : lights) {
          color.red = color.red + ambientReflection.red * i.ambientLight.red / 255;
          color.blue = color.blue + ambientReflection.blue * i.ambientLight.blue /255;
          color.green = color.green + ambientReflection.green * i.ambientLight.green /255;

        }
        image(xi, yi) = color;
        buffer(xi, yi) = oz;
      }
    }
  }
  

}


double figures_3d::min(std::vector<double> values) {
  double min = values[0];
  for (auto i : values){
    if (min > i ){
      min = i;
    }
  }
  return min;
}

double figures_3d::max(std::vector<double> values) {
  double max = values[0];
  for (auto i : values){
    if (max < i ){
      max = i;
    }
  }
  return max;
}

bool in_near_far() {

}

figures_3d::Figure figures_3d::clipping(Figure fig, const ini::Configuration &conf) {

  double near = conf["General"]["dNear"].as_double_or_die();
  double far = conf["General"]["dFar"].as_double_or_die();
  double hfov = conf["General"]["hfov"].as_double_or_die() * M_PI /180;
  double ratio = conf["General"]["aspectRatio"].as_double_or_die();
  double right = near * tan(hfov/2);
  double top = right/ ratio;
  ini::DoubleTuple dir = conf["General"]["viewDirection"].as_double_tuple_or_die();

  fig.faces = triangulate(fig.faces);
  //todo near clipping
  
  for (auto j = 0; j < fig.faces.size(); j++) {
    std::vector<int> in;
    std::vector<int> out;
    //
    for (auto &k : fig.faces[j]) {
      if (fig.points[k].z <= - near ) {
        in.emplace_back(k);
      } else {
        out.emplace_back(k);
      }
    }
    //
    if (in.size() == 0) {
      fig.faces.erase(fig.faces.begin() + j);
      j--;
    } else if (in.size() == 2) {
      auto & A = fig.points[out[0]];
      auto & B = fig.points[in[0]];
      auto & C = fig.points[in[1]];

      double p1 = (-near - B.z) / (A.z - B.z);
      double p2 = (-near - C.z) / (A.z - C.z);

      fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
      fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

      fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 2), (int) (fig.points.size() - 1), in[1]});
      fig.faces.erase(fig.faces.begin() + j);
      j--;

    } else if (in.size() == 1) {

      auto & A = fig.points[in[0]];
      auto & B = fig.points[out[0]];
      auto & C = fig.points[out[1]];

      double p1 = (-near - B.z) / (A.z - B.z);
      double p2 = (-near - C.z) / (A.z - C.z);

      fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
      fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

      fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 1), (int) (fig.points.size() - 2)});
      fig.faces.erase(fig.faces.begin() + j);
      j--;
    } else { //doe niks
    }
  }
    fig.faces = triangulate(fig.faces);

    //todo far clipping
    for (auto j = 0; j < fig.faces.size(); j++) {
      std::vector<int> in;
      std::vector<int> out;
      //
      for (auto &k : fig.faces[j]) {
        if (fig.points[k].z >= - far ) {
          in.emplace_back(k);
        } else {
          out.emplace_back(k);
        }
      }
      //
      if (in.size() == 0) {
        fig.faces.erase(fig.faces.begin() + j);
        j--;
      } else if (in.size() == 2) {
        auto & A = fig.points[out[0]];
        auto & B = fig.points[in[0]];
        auto & C = fig.points[in[1]];

        double p1 = (far - B.z) / (A.z - B.z);
        double p2 = (far - C.z) / (A.z - C.z);

        fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
        fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

        fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 2), (int) (fig.points.size() - 1), in[1]});
        fig.faces.erase(fig.faces.begin() + j);
        j--;

      } else if (in.size() == 1) {

        auto & A = fig.points[in[0]];
        auto & B = fig.points[out[0]];
        auto & C = fig.points[out[1]];

        double p1 = (far - B.z) / (A.z - B.z);
        double p2 = (far - C.z) / (A.z - C.z);

        fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
        fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

        fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 1), (int) (fig.points.size() - 2)});
        fig.faces.erase(fig.faces.begin() + j);
        j--;
      } else { //doe niks
      }
    }
    fig.faces = triangulate(fig.faces);

    //todo right clipping
    int faces_size_now = fig.faces.size();
    for (auto j = 0; j < fig.faces.size(); j++) {
      std::vector<int> in;
      std::vector<int> out;
      //
      for (auto &k : fig.faces[j]) {
        if (fig.points[k].x * near / -fig.points[k].z <= right * 1.002) {
          in.emplace_back(k);
        } else {
          out.emplace_back(k);
        }
      }
      //
      if (in.size() == 0) {
        fig.faces.erase(fig.faces.begin() + j);
        j--;
      } else if (in.size() == 2) {
        
        auto & A = fig.points[out[0]];
        auto & B = fig.points[in[0]];
        auto & C = fig.points[in[1]];
        
        double p1 = (B.x * near + B.z * right) / ((B.x - A.x) * near + (B.z - A.z) * right);
        double p2 = (C.x * near + C.z * right) / ((C.x - A.x) * near + (C.z - A.z) * right);

        fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
        fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

        fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 2), (int) (fig.points.size() - 1), in[1]});
        fig.faces.erase(fig.faces.begin() + j);
        j--;

      } else if (in.size() == 1) {
        
        auto & A = fig.points[in[0]];
        auto & B = fig.points[out[0]];
        auto & C = fig.points[out[1]];

        double p1 = (B.x * near + B.z * right) / ((B.x - A.x) * near + (B.z - A.z) * right);
        double p2 = (C.x * near + C.z * right) / ((C.x - A.x) * near + (C.z - A.z) * right);

        fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
        fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

        fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 1), (int) (fig.points.size() - 2)});
        fig.faces.erase(fig.faces.begin() + j);
        j--;
      } else { // niks doen
      }
    }
    fig.faces = triangulate(fig.faces);

    //todo  left clipping
    for (auto j = 0; j < fig.faces.size(); j++) {
      std::vector<int> in;
      std::vector<int> out;
      //
      for (auto &k : fig.faces[j]) {
        if (fig.points[k].x * near / -fig.points[k].z >= -right * 1.002) {
          in.emplace_back(k);
        } else {
          out.emplace_back(k);
        }
      }
      //
      if (in.size() == 0) {
        fig.faces.erase(fig.faces.begin() + j);
        j--;
      } else if (in.size() == 2) {

        auto & A = fig.points[out[0]];
        auto & B = fig.points[in[0]];
        auto & C = fig.points[in[1]];

        double p1 = (B.x * near + B.z * -right) / ((B.x - A.x) * near + (B.z - A.z) * -right);
        double p2 = (C.x * near + C.z * -right) / ((C.x - A.x) * near + (C.z - A.z) * -right);

        fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
        fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

        fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 2), (int) (fig.points.size() - 1), in[1]});
        fig.faces.erase(fig.faces.begin() + j);
        j--;

      } else if (in.size() == 1) {

        auto & A = fig.points[in[0]];
        auto & B = fig.points[out[0]];
        auto & C = fig.points[out[1]];

        double p1 = (B.x * near + B.z * -right) / ((B.x - A.x) * near + (B.z - A.z) * -right);
        double p2 = (C.x * near + C.z * -right) / ((C.x - A.x) * near + (C.z - A.z) * -right);

        fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
        fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

        fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 1), (int) (fig.points.size() - 2)});
        fig.faces.erase(fig.faces.begin() + j);
        j--;
      } else { // niks doen
      }
    }
    fig.faces = triangulate(fig.faces);


    //todo  top clipping
    for (auto j = 0; j < fig.faces.size(); j++) {
      std::vector<int> in;
      std::vector<int> out;
      //
      for (auto &k : fig.faces[j]) {
        if (fig.points[k].y * near / -fig.points[k].z <= top * 1.002 ) {
          in.emplace_back(k);
        } else {
          out.emplace_back(k);
        }
      }
      //
      if (in.size() == 0) {
        fig.faces.erase(fig.faces.begin() + j);
        j--;
      } else if (in.size() == 2) {

        auto & A = fig.points[out[0]];
        auto & B = fig.points[in[0]];
        auto & C = fig.points[in[1]];

        double p1 = (B.y * near + B.z * top) / ((B.y - A.y) * near + (B.z - A.z) * top);
        double p2 = (C.y * near + C.z * top) / ((C.y - A.y) * near + (C.z - A.z) * top);

        fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
        fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

        fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 2), (int) (fig.points.size() - 1), in[1]});
        fig.faces.erase(fig.faces.begin() + j);
        j--;

      } else if (in.size() == 1) {

        auto & A = fig.points[in[0]];
        auto & B = fig.points[out[0]];
        auto & C = fig.points[out[1]];

        double p1 = (B.y * near + B.z * top) / ((B.y - A.y) * near + (B.z - A.z) * top);
        double p2 = (C.y * near + C.z * top) / ((C.y - A.y) * near + (C.z - A.z) * top);

        fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
        fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

        fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 1), (int) (fig.points.size() - 2)});
        fig.faces.erase(fig.faces.begin() + j);
        j--;
      } else { // niks doen
      }
    }
    fig.faces = triangulate(fig.faces);

    //todo bottom clipping

    for (auto j = 0; j < fig.faces.size(); j++) {
      std::vector<int> in;
      std::vector<int> out;
      //
      for (auto &k : fig.faces[j]) {
        if (fig.points[k].y * near / -fig.points[k].z >= - top * 1.002 ) {
          in.emplace_back(k);
        } else {
          out.emplace_back(k);
        }
      }
      //
      if (in.size() == 0) {
        fig.faces.erase(fig.faces.begin() + j);
        j--;
      } else if (in.size() == 2) {

        auto & A = fig.points[out[0]];
        auto & B = fig.points[in[0]];
        auto & C = fig.points[in[1]];

        double p1 = (B.y * near + B.z * -top) / ((B.y - A.y) * near + (B.z - A.z) * -top);
        double p2 = (C.y * near + C.z * -top) / ((C.y - A.y) * near + (C.z - A.z) * -top);

        fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
        fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

        fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 2), (int) (fig.points.size() - 1), in[1]});
        fig.faces.erase(fig.faces.begin() + j);
        j--;

      } else if (in.size() == 1) {

        auto & A = fig.points[in[0]];
        auto & B = fig.points[out[0]];
        auto & C = fig.points[out[1]];

        double p1 = (B.y * near + B.z * -top) / ((B.y - A.y) * near + (B.z - A.z) * -top);
        double p2 = (C.y * near + C.z * -top) / ((C.y - A.y) * near + (C.z - A.z) * -top);

        fig.points.emplace_back(Vector3D::point(p1 * A + (1-p1) * B));
        fig.points.emplace_back(Vector3D::point(p2 * A + (1-p2) * C));

        fig.faces.emplace_back(Face{in[0], (int) (fig.points.size() - 1), (int) (fig.points.size() - 2)});
        fig.faces.erase(fig.faces.begin() + j);
        j--;
      } else { // niks doen
      }
    }
    fig.faces = triangulate(fig.faces);

  return fig;
}

figures_3d::Figures3D::Figures3D(figures_3d::Figure fig) {
  this->push_back(fig);
}
