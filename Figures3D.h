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
#include <fstream>

#include "vector3d.h"
#include "easy_image.h"
#include "lines.h"
#include "ZBuffer.h"
#include "l_system.h"
#include "Light.h"


namespace figures_3d {

    typedef std::vector<int> Face;

    class Figure {
     public:

        std::vector<Vector3D> points;
        std::vector<Face> faces;
        img::Color ambientReflection;
        img::Color diffuseReflection;
        img::Color specularReflection;
        double reflectionCoefficient;

        Figure() = default;

        Figure(std::vector<Vector3D> points, std::vector<Face> faces, img::Color &color)
                    : points(std::move(std::move(points))), faces(std::move(std::move(faces))), ambientReflection(color) {}

        void scale(const double &scale);

        void rotateX(const double &angle);

        void rotateY(const double &angle);

        void rotateZ(const double &angle);

        void translate(const double &, const double &, const double &);

        void to_eye( const double &theta, const double &phi, const double &r);

        void to_eye_clip(const double &theta, const double &phi, Vector3D &);

        lines_2d::Lines2D project_figure(double d = 1);

        img::Color calcColor(Lights3D &);


    };

    //typedef std::vector<Figure> Figures3D;

    struct Figures3D :  public std::vector<Figure> {

      Figures3D(figures_3d::Figure fig);

      Figures3D() = default;

    };

    figures_3d::Figure clipping(Figure figs, const ini::Configuration &conf);

    figures_3d::Figure calc_fig(LParser::LSystem3D&);

    figures_3d::Figures3D fractal(int ,const double&, Figure);

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

    void draw_triangle_light(const Vector3D &A,
                             const Vector3D &B,
                             const Vector3D &C,
                             double d,
                             double dx,
                             double dy,
                             img::EasyImage &image,
                             zbuffer::ZBuffer &buffer,
                             img::Color ambientReflection,
                             img::Color diffuseReflection,
                             img::Color specularReflection, double reflectionCoeff,
                             Lights3D& lights);

}




#endif //ENGINE_FIGURES3D_H
