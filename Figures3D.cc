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
        fig.points.push_back(Vector3D::point(1, 1, 1));
        fig.points.push_back(Vector3D::point(1, 1, -1));
        fig.points.push_back(Vector3D::point(1, -1, 1));
        fig.points.push_back(Vector3D::point(1, -1, -1));
        fig.points.push_back(Vector3D::point(-1, 1, 1));
        fig.points.push_back(Vector3D::point(-1, 1, -1));
        fig.points.push_back(Vector3D::point(-1, -1, 1));
        fig.points.push_back(Vector3D::point(-1, -1, -1));
        fig.faces.push_back(std::vector<int>{0, 1});
        fig.faces.push_back(std::vector<int>{0, 4});
        fig.faces.push_back(std::vector<int>{0, 2});
        fig.faces.push_back(std::vector<int>{2, 6});
        fig.faces.push_back(std::vector<int>{2, 3});
        fig.faces.push_back(std::vector<int>{4, 6});
        fig.faces.push_back(std::vector<int>{4, 5});
        fig.faces.push_back(std::vector<int>{6, 7});
        fig.faces.push_back(std::vector<int>{1, 3});
        fig.faces.push_back(std::vector<int>{1, 5});
        fig.faces.push_back(std::vector<int>{3, 7});
        fig.faces.push_back(std::vector<int>{5, 7});
        return fig;
    }

    inline figures_3d::Figure cylinder(const ini::Configuration &configuration, unsigned int i) {
        figures_3d::Figure fig;
        double height = configuration["Figure" + std::to_string(i)]["height"].as_double_or_die();
        int points_size = configuration["Figure" + std::to_string(i)]["n"].as_int_or_die();

        auto first = Vector3D::point(1, 0, 0); fig.points.push_back(first);
        auto second = Vector3D::point(1, 0, height); fig.points.push_back(second);
        std::vector<int> firts_face = {0, 1};

        for (auto j = 2; j < points_size *2 +2 ; j++){ //todo +2 lost iets op
            auto point1 = Vector3D::point(std::cos(2 * M_PI * j  / points_size), std::sin(2 * M_PI * j  / points_size), 0);
            auto point2 = Vector3D::point(std::cos(2 * M_PI * j  / points_size), std::sin(2 * M_PI * j  / points_size), height);
            fig.points.insert(fig.points.begin(), point1);
            fig.points.insert(fig.points.begin(), point2);
            std::vector<int> face1;
            face1.push_back(j); face1.push_back(j + 2);
            fig.faces.push_back(face1);
            if (j % 2 == 1){
                std::vector<int> face2; face2.push_back(j); face2.push_back(j - 1);
                fig.faces.push_back(face2);
            }
        }

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
            for (auto j = 0; j < m ; j++){ // laatste niet meetellen , gaat dan scheef , dus heeft geen zin
                std::vector<int> face1 = {i * n + j, i * n + j + 1};
                std::vector<int> face2 = {i * n + j, i * n + j + m};

                /// {vreemde shit} maar het werkt
                if (j == m-1){
                    std::vector<int> face3 = {i * n + j, i * n};
                    fig.faces.push_back(face3);
                }
                else{
                    fig.faces.push_back(face1);
                }
                if (i != n-1){
                    fig.faces.push_back(face2);
                }
                else{
                    std::vector<int> face4 = {i * n + j, j};
                    fig.faces.push_back(face4);
                }
                // {vreemde shit}
            }
//            std::vector<int> face_last = {i * n , i * n + m-1};
            //fig.faces.push_back(face_last);
        }
        return fig;
    }

    inline figures_3d::Figure tetrahedron(const ini::Configuration &configuration, unsigned int i){
        figures_3d::Figure fig;
        fig.points.push_back(Vector3D::point(1, -1, -1));
        fig.points.push_back(Vector3D::point(-1, 1, -1));
        fig.points.push_back(Vector3D::point(1, 1, 1));
        fig.points.push_back(Vector3D::point(-1, -1, 1));
        for (auto i = 0; i < 4; i++ ){
            for (auto j = i + 1 ; j < 4; j++ ){
                std::vector<int> face = { i, j};
                fig.faces.push_back(face);
            }
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
        for (auto i = 0; i < this->faces[k].size(); i ++){
            for (auto j = i+1; j < this->faces[k].size() ; j++ ){
                Vector3D &point1 = this->points[this->faces[k][i]];
                Vector3D &point2 = this->points[this->faces[k][j]];
                lines.push_back(lines_2d::Line2D(
                        lines_2d::Point2D(- point1.x * d / point1.z , - point1.y * d / point1.z),
                        lines_2d::Point2D(- point2.x * d / point2.z , - point2.y * d / point2.z),
                        this->color ));
            }
        }
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


        ini:: DoubleTuple fig_color = configuration["Figure" + std::to_string(i)]["color"].as_double_tuple_or_die();
        figure.color = img::Color(fig_color[0]*255, fig_color[1]*255, fig_color[2]*255);
        double r = std::pow(std::pow(eye[0], 2.0) + std::pow(eye[1], 2.0) + std::pow(eye[2], 2.0), 0.5);
        ini::DoubleTuple center = configuration["Figure" + std::to_string(i)]["center"].as_double_tuple_or_die();
        figure.translate(center[0], center[1], center[2]);
        figure.scale(configuration["Figure" + std::to_string(i)]["scale"].as_double_or_die());
        figure.rotateX(M_PI /180 * configuration["Figure" + std::to_string(i)]["rotateX"].as_int_or_die());
        figure.rotateY(M_PI /180 * configuration["Figure" + std::to_string(i)]["rotateY"].as_int_or_die());
        figure.rotateZ(M_PI /180 * configuration["Figure" + std::to_string(i)]["rotateZ"].as_int_or_die());
        figure.to_eye( std::atan2(eye[1],eye[0]), std::acos(eye[2] / r), r);
        figs.push_back(figure);
    }
    return figs;
}
