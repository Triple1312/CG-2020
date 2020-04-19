//
// Created by Phili on 3/4/2020.
//

#ifndef ENGINE_FIGURES3D_H
#define ENGINE_FIGURES3D_H

#include <iostream>
#include <utility>
#include <utility>
#include <vector>
#include <list>
#include <cmath>
#include <tuple>

#include "vector3d.h"
#include "easy_image.h"
#include "lines.h"
#include "ZBuffer.h"

namespace figures_3d {

    typedef std::vector<int> Face;

    struct Figure {

        std::vector<Vector3D> points;
        std::vector<Face> faces;
        img::Color color;

        Figure() = default;

        Figure(std::vector<Vector3D> points, std::vector<Face> faces, img::Color &color)
                    : points(std::move(std::move(points))), faces(std::move(std::move(faces))), color(color) {}

        void scale(const double &scale);

        void rotateX(const double &angle);

        void rotateY(const double &angle);

        void rotateZ(const double &angle);

        void translate(const double &, const double &, const double &);

        void to_eye( const double &theta, const double &phi, const double &r);

        lines_2d::Lines2D project_figure(double d = 1);
    };

    typedef std::vector<Figure> Figures3D;

    Figures3D Wireframe(const ini::Configuration &config);

    std::vector<figures_3d::Face> triangulate(Face&);

    std::vector<figures_3d::Face> triangulate(std::vector<figures_3d::Face>&);

    double min(std::vector<double>);

    double max(std::vector<double>);

    void draw_triangle(const Vector3D &A,
                       const Vector3D &B,
                       const Vector3D &C,
                       double dx,
                       double dy,
                       img::EasyImage &image,
                       zbuffer::ZBuffer &buffer,
                       img::Color &color,
                       double d = 1);



}



#endif //ENGINE_FIGURES3D_H
