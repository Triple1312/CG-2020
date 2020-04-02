//
// Created by Phili on 2/27/2020.
//

#ifndef LINES_H_
#define LINES_H_

#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include <utility>
#include <tuple>

#include "l_parser.h"
#include "ini_configuration.h"
#include "easy_image.h"


/// a namespace for functions of 2d lines
namespace lines_2d {

    /// struct that defines a point of a line
    struct Point2D {
        Point2D(double x, double y) : x(x), y(y) {}

        explicit Point2D(std::pair<double, double> coords)
            : x(coords.first), y(coords.second) {}
        double x;
        double y;

        Point2D();
    };

    /// struct that defines a line
    struct Line2D {
        Line2D(Point2D p1, Point2D p2, img::Color c)
            : p1(p1), p2(p2), color(c) {}
        Point2D p1;
        Point2D p2;
        img::Color color;

        Line2D();
    };

    using Lines2D = std::vector<Line2D>;


    /// Function that returns a pair with
    /// the minimum x and maximum x of a vector of lines.
    /// \param lines
    /// \return minimum and maximum of a vector of lines.
    std::pair<double, double> min_max_x(const Lines2D& lines);

    /// Function that returns a pair with
    /// the minimum x and maximum x of a vector of lines.
    /// \param lines
    /// \return
    std::pair<double, double> min_max_y(const Lines2D& lines);

    /// Function that scales all the lines to fit in the given size.
    /// \param lines
    /// \param size
    /// \return std::pair<double, double> with the size for x and y of picture
    std::pair<int, int> scaleLines(
            Lines2D &lines, const int& size);

}  // namespace lines_2d

namespace lines_3d{

    struct Vector3D{
        double x;
        double y;
        double z;

        Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}
    };

    struct Point3D{
        double x;
        double y;
        double z;

        Point3D(double x, double y, double z) : x(x), y(y), z(z) {}

        explicit Point3D(std::vector<double>);  //kan effectiever

        explicit Point3D(std::tuple<double, double, double> coords);

        Vector3D operator - (const Point3D& point);

    };






}

#endif  // LINES_H_