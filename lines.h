//
// Created by Phili on 2/27/2020.
//

#ifndef ENGINE_LINES_H
#define ENGINE_LINES_H

#include <iostream>
#include <string>
#include <memory>

#include "l_parser.h"
#include "ini_configuration.h"
#include "easy_image.h"

namespace lines_2d{

    struct Point2D{

        Point2D(double x, double y) : x(x), y(y){};

        explicit Point2D(std::pair<double, double> coords) : x(coords.first), y(coords.second){};

        double x;
        double y;
    };

    struct Line2D{

        Line2D(Point2D p1, Point2D p2, img::Color c) : p1(p1), p2(p2), color(c){};

        Point2D p1;
        Point2D p2;
        img::Color color;
    };

    using Lines2D = std::vector<Line2D>;

    std::shared_ptr<img::EasyImage> draw2DLines(Lines2D& lines, int size);

    /// Function that returns a pair with the minimum x and maximum x of a vector of lines.
    /// \param lines
    /// \return
    std::pair<double, double> min_max_x(Lines2D& lines);

    /// Function that returns a pair with the minimum x and maximum x of a vector of lines.
    /// \param lines
    /// \return
    std::pair<double, double> min_max_y(Lines2D& lines);

    /// Function that scales all the lines to fit in the given size.
    /// \param lines
    /// \param size
    /// \return std::pair<double, double> with the size for x and y of picture
    std::pair<double, double> scaleLines(Lines2D& lines, const int& size);

    ///
    /// \param lines
    /// \return the ange of the x coordinates of lines
    double rangeX(Lines2D& lines);

    ///
    /// \param lines
    /// \return the ange of the y coordinates of lines
    double rangeY(Lines2D& lines);

} // namespace lines2D



#endif //ENGINE_LINES_H
