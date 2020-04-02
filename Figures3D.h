//
// Created by Phili on 3/4/2020.
//

#ifndef ENGINE_FIGURES3D_H
#define ENGINE_FIGURES3D_H

#include <iostream>
#include <vector>
#include <list>
#include <cmath>

#include "vector3d.h"
#include "easy_image.h"
#include "lines.h"

namespace figures_3d {

    typedef std::vector<int> Face;

    struct Figure {

        std::vector<Vector3D> points;
        std::vector<Face> faces;
        img::Color color;

        Figure() {}

        Figure(std::vector<Vector3D> points, std::vector<Face> faces, img::Color &color)
                    : points(points), faces(faces), color(color) {}

        void scale(const double &scale);

        void rotateX(const double &angle);

        void rotateY(const double &angle);

        void rotateZ(const double &angle);

        void translate(const double &, const double &, const double &);

        void to_eye( const double &theta, const double &phi, const double &r);

        lines_2d::Lines2D project_figure(double d = 1);
    };

    typedef std::vector<Figure> Figures3D;






}



#endif //ENGINE_FIGURES3D_H
