//
// Created by Phili on 3/4/2020.
//

#include "Figures3D.h"

namespace {

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
    scaler( 0, 0 ) = scale;
    scaler( 1, 1 ) = scale;
    scaler( 2, 2 ) = scale;
    scaler( 3, 3 ) = 1;

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


