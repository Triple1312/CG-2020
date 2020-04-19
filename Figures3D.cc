//
// Created by Phili on 3/4/2020.
//

#include "Figures3D.h"

namespace {

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

    inline figures_3d::Figure cube(const ini::Configuration &configuration, unsigned int i){
        figures_3d::Figure fig;
        fig.points.push_back(Vector3D::point(1, -1, -1));
        fig.points.push_back(Vector3D::point(-1, 1, -1));
        fig.points.push_back(Vector3D::point(1, 1, 1));
        fig.points.push_back(Vector3D::point(-1, -1, 1));
        fig.points.push_back(Vector3D::point(1, 1, -1));
        fig.points.push_back(Vector3D::point(-1, -1, -1));
        fig.points.push_back(Vector3D::point(1, -1, 1));
        fig.points.push_back(Vector3D::point(-1, 1, 1));
        fig.faces.push_back(std::vector<int>{0, 5, 1, 4});
        fig.faces.push_back(std::vector<int>{6, 3, 7, 2});
        fig.faces.push_back(std::vector<int>{1, 5, 3, 7});
        fig.faces.push_back(std::vector<int>{2, 6, 0, 4});
        fig.faces.push_back(std::vector<int>{0, 5, 3, 6});
        fig.faces.push_back(std::vector<int>{2, 7, 1, 4});
        return fig;
    }

    inline figures_3d::Figure cylinder(const ini::Configuration &configuration, unsigned int i) {
        figures_3d::Figure fig;
        double height = configuration["Figure" + std::to_string(i)]["height"].as_double_or_die();
        int n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
        std::vector<int> firts_face = {0, 1};

        for (int i = 0; i < n; i++) {
          fig.points.emplace_back(Vector3D::point(std::cos(2 * M_PI * i  / n),
                                                  std::sin(2 * M_PI * i  / n), 0));
        }
        for (int i = 0; i < n; i++) {
          fig.points.emplace_back(Vector3D::point(std::cos(2 * M_PI * i  / n),
                                                  std::sin(2 * M_PI * i  / n), height));
        }
        for (int j = 0; j < n-1; j++) {
          fig.faces.push_back(std::vector<int>{j, j+1, j+n+1, j+n});
        }
        fig.faces.push_back(std::vector<int>{0, n-1, 2*n -1, n});
        std::vector<int> onder;
        std::vector<int> boven;
        for (int k = 0; k < n; k++) {
          onder.push_back(k);
          boven.push_back(n+k);
        }
        fig.faces.push_back(onder);
        fig.faces.push_back(boven);
        return fig;
    }

    inline figures_3d::Figure cone(const ini::Configuration &configuration, unsigned int i) {
        figures_3d::Figure fig;
        double height = configuration["Figure" + std::to_string(i)]["height"].as_double_or_die();
        int points_size = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
        auto top = Vector3D::point(0, 0, height);
        fig.points.push_back(top);
        auto start_point = Vector3D::point(1, 0, 0); fig.points.push_back(start_point);

        for (auto j = 1; j <= points_size; j++){
            auto point = Vector3D::point(std::cos(2 * M_PI * j  / points_size), std::sin(2 * M_PI * j  / points_size), 0);
            fig.points.push_back(point);
            std::vector<int> face;
             face.push_back(j); face.push_back(j - 1); face.push_back(0);
            fig.faces.push_back(face);
        }
        std::vector<int> last_face = {0, 1, int(fig.points.size()-2)};
        fig.faces.push_back(last_face);
        return fig;
    }

    inline figures_3d::Figure torus(const ini::Configuration &configuration, unsigned int i){
        figures_3d::Figure fig;
        auto R = configuration["Figure" + std::to_string(i)]["R"].as_double_or_die();
        auto r = configuration["Figure" + std::to_string(i)]["r"].as_double_or_die();
        auto n = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();
        auto m = configuration["Figure" + std::to_string(i)]["m"].as_int_or_die();

        for (auto i = 0; i < n; i++){
            for (auto j = 0; j < m; j++){
                double u = 2 * i * M_PI / n;
                double v = 2 * j * M_PI / m;
                fig.points.push_back(Vector3D::point(((R + r * std::cos(v)) * std::cos(u)),
                                                     ((R + r * std::cos(v)) * std::sin(u)),
                                                     r * std::sin(v)));
            }
        }
        for (auto i = 0; i < n; i++){
          for (auto j = 0; j < m; j++){
            fig.faces.push_back(std::vector<int>{i*m+j, ((i+1)%n)*m+j,((i+1)%n)*m+((j+1)%m), i*m+((j+1)%m)});
          }
        }
        return fig;
    }

    inline figures_3d::Figure tetrahedron(const ini::Configuration &configuration, unsigned int i){
        figures_3d::Figure fig;
        fig.points.push_back(Vector3D::point(1, -1, -1));
        fig.points.push_back(Vector3D::point(-1, 1, -1));
        fig.points.push_back(Vector3D::point(1, 1, 1));
        fig.points.push_back(Vector3D::point(-1, -1, 1));
        fig.faces.push_back(std::vector<int>{0, 1, 2});
        fig.faces.push_back(std::vector<int>{0, 2, 3});
        fig.faces.push_back(std::vector<int>{0, 1, 3});
        fig.faces.push_back(std::vector<int>{1, 2, 3});
        return fig;
    }

    inline std::vector<std::vector<int>> ico_faces(){
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

    inline figures_3d::Figure icosahedron() {
        figures_3d::Figure fig;
        fig.points.push_back(Vector3D::point(0, 0, sqrt(5) / 2));
        for (auto i = 2; i <= 6; i++) {
            fig.points.push_back(Vector3D::point(std::cos((i-2) * 2 * M_PI / 5),
                                                 std::sin((i-2) * 2 * M_PI / 5),
                                                 0.5));
        }
        for (int j = 7; j <= 11 ; j++) {
            fig.points.push_back(Vector3D::point(std::cos( M_PI / 5 + (j-7) * 2 * M_PI / 5),
                                                 std::sin( M_PI / 5 + (j-7) * 2 * M_PI / 5),
                                                 -0.5));
        }
        fig.points.push_back(Vector3D::point(0, 0, - sqrt(5) / 2));

        fig.faces = ico_faces();
        return fig;
    }

    inline figures_3d::Figure octahedron() {
        figures_3d::Figure fig;
        fig.points.push_back(Vector3D::point(1, 0, 0 ));
        fig.points.push_back(Vector3D::point(0, 1, 0 ));
        fig.points.push_back(Vector3D::point(-1, 0, 0 ));
        fig.points.push_back(Vector3D::point(0, -1, 0 ));
        fig.points.push_back(Vector3D::point(0, 0, -1 ));
        fig.points.push_back(Vector3D::point(0, 0, 1 ));

        fig.faces.push_back(std::vector<int>{1, 2, 6});
        fig.faces.push_back(std::vector<int>{2, 3, 6});
        fig.faces.push_back(std::vector<int>{3, 4, 6});
        fig.faces.push_back(std::vector<int>{4, 1, 6});
        fig.faces.push_back(std::vector<int>{2, 1, 5});
        fig.faces.push_back(std::vector<int>{3, 2, 5});
        fig.faces.push_back(std::vector<int>{4, 3, 5});
        fig.faces.push_back(std::vector<int>{1, 4, 5});

        return fig;
    }

    inline figures_3d::Figure dodecahedron(){
        figures_3d::Figure ico = icosahedron();
        figures_3d::Figure fig;

        for (auto i : ico.faces){
            auto &points = ico.points;
            fig.points.push_back(Vector3D::point(( points[i[0]].x + points[i[1]].x + points[i[2]].x) / 3,
                                                 ( points[i[0]].y + points[i[1]].y + points[i[2]].y) / 3,
                                                 ( points[i[0]].z + points[i[1]].z + points[i[2]].z) / 3) );
        }

        for (int j = 0; j < fig.points.size(); j++) {
            for (int k = j+1; k < fig.points.size(); k++) {
                if ( sqrt(pow(fig.points[j].x + fig.points[k].x, 2) +
                          pow(fig.points[j].y + fig.points[k].y, 2) +
                          pow(fig.points[j].z + fig.points[k].z, 2)) > 1.5){
                    fig.faces.push_back(std::vector<int>{j, k}); // todo fixen
                }
            }
        }

        return fig;
    }

    inline figures_3d::Figure sphere(int n){
        figures_3d::Figure fig = icosahedron();

        for (auto j = 0; j < n; j++) {
            auto faces_b_s = fig.faces.size(); // begin size of faces
            for (int i = 0; i < faces_b_s; i++) {
                fig.points.push_back((fig.points[fig.faces[0][0]] + fig.points[fig.faces[0][1]]) / 2);
                fig.points.push_back((fig.points[fig.faces[0][0]] + fig.points[fig.faces[0][2]]) / 2);
                fig.points.push_back((fig.points[fig.faces[0][1]] + fig.points[fig.faces[0][2]]) / 2);

                fig.faces.push_back(
                        std::vector<int>{fig.faces[0][0], int(fig.points.size()) - 3, int(fig.points.size()) - 2});
                fig.faces.push_back(
                        std::vector<int>{fig.faces[0][1], int(fig.points.size()) - 3, int(fig.points.size()) - 1});
                fig.faces.push_back(
                        std::vector<int>{fig.faces[0][2], int(fig.points.size()) - 2, int(fig.points.size()) - 1});
                fig.faces.push_back(std::vector<int>{int(fig.points.size()) - 1, int(fig.points.size()) - 3,
                                                     int(fig.points.size()) - 2});
                fig.faces.erase(fig.faces.begin());
            }

        }

        for (int k = 0; k < fig.points.size(); k++) {
            double r = sqrt(pow(fig.points[k].x, 2) +
                            pow(fig.points[k].y, 2) +
                            pow(fig.points[k].z, 2));
            fig.points[k] = fig.points[k] / r;
        }
        return fig;
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

    inline Matrix getMatrix_translate(const double &x, const double &y, const double &z){
        Matrix translate_matrix;
        translate_matrix( 4, 1) = x;
        translate_matrix( 4, 2) = y;
        translate_matrix( 4, 3) = z;
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

//    inline lines_2d::Point2D point3D_to_2D(const Vector3D &point){
//        return lines_2d::Point2D(std::)
//    }
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
                    this->color, point1.z, point2.z ));

        }
        Vector3D &point1 = this->points[this->faces[k][0]];
        Vector3D &point2 = this->points[this->faces[k][this->faces[k].size()-1]];
        lines.push_back(lines_2d::Line2D(
            lines_2d::Point2D(- point1.x * d / point1.z , - point1.y * d / point1.z),
            lines_2d::Point2D(- point2.x * d / point2.z , - point2.y * d / point2.z),
            this->color, point1.z, point2.z ));
    }
    return lines;
}

void figures_3d::Figure::to_eye( const double &theta, const double &phi, const double &r) {
    for (auto &i : this->points){
        i *= eye_matrix(theta, phi, r);
    }
}


figures_3d::Figures3D figures_3d::Wireframe(const ini::Configuration &configuration) {
    figures_3d::Figures3D figs;
    int blub = configuration["General"]["nrFigures"].as_int_or_die();

    ini::DoubleTuple eye = configuration["General"]["eye"].as_double_tuple_or_die();
    for (auto i = 0; i < configuration["General"]["nrFigures"].as_int_or_die(); i++){

        figures_3d::Figure figure;  // todo kan veel beter door direct to initializen in de vector
        if (configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "LineDrawing"){
            figure = line_drawing(configuration, i);  // i is the number of the figure
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Cube"){
            figure = cube(configuration, i);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Cone"){
            figure = cone(configuration, i);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Cylinder"){
            figure = cylinder(configuration, i);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Torus"){
            figure = torus(configuration, i);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Tetrahedron"){
            figure = tetrahedron(configuration, i);
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Icosahedron"){
            figure = icosahedron();
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Octahedron"){
            figure = octahedron();
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Dodecahedron"){
            figure = dodecahedron();
        }
        else if(configuration["Figure" + std::to_string(i)]["type"].as_string_or_die() == "Sphere"){
            figure = sphere(configuration["Figure" + std::to_string(i)]["n"].as_int_or_die());
        }

        ini:: DoubleTuple fig_color = configuration["Figure" + std::to_string(i)]["color"].as_double_tuple_or_die();
        figure.color = img::Color(fig_color[0]*255, fig_color[1]*255, fig_color[2]*255);
        double r = std::pow(std::pow(eye[0], 2.0) + std::pow(eye[1], 2.0) + std::pow(eye[2], 2.0), 0.5);
        ini::DoubleTuple center = configuration["Figure" + std::to_string(i)]["center"].as_double_tuple_or_die();

        figure.scale(configuration["Figure" + std::to_string(i)]["scale"].as_double_or_die());

        figure.rotateX(M_PI /180 * configuration["Figure" + std::to_string(i)]["rotateX"].as_double_or_die());
        figure.rotateY(M_PI /180 * configuration["Figure" + std::to_string(i)]["rotateY"].as_double_or_die());
        figure.rotateZ(M_PI /180 * configuration["Figure" + std::to_string(i)]["rotateZ"].as_double_or_die());
        figure.translate(center[0], center[1], center[2]);
        figure.to_eye( std::atan2(eye[1],eye[0]), std::acos(eye[2] / r), r);
        figs.push_back(figure);
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
