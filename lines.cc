//
// Created by Phili on 2/27/2020.
//

#include <cmath>
#include "lines.h"

namespace lines_2d{

    std::pair<double, double> min_max_y(Lines2D& lines){
        double ymin = lines[0].p1.y;
        double ymax = lines[0].p1.y;
        for (auto yp : lines){
            if (yp.p1.y < ymin){ ymin = yp.p1.y;}
            if (yp.p1.y > ymax){ ymax = yp.p1.y;}
            if (yp.p2.y < ymin){ ymin = yp.p2.y;}
            if (yp.p2.y > ymax){ ymax = yp.p2.y;}
        }
        return std::pair<double,double>(ymin,ymax);
    }

    std::pair<double, double> min_max_x(Lines2D& lines){
        double xmin = lines[0].p1.x;
        double xmax = lines[0].p1.x;
        for (auto xp : lines){  // kan verbeterd worden
            if (xp.p1.x < xmin){ xmin = xp.p1.x;}
            if (xp.p1.x > xmax){ xmax = xp.p1.x;}
            if (xp.p2.x < xmin){ xmin = xp.p2.x;}
            if (xp.p2.x > xmax){ xmax = xp.p2.x;}
        }
        return std::pair<double,double>(xmin,xmax);
    }

    double rangeX(Lines2D& lines){
        auto minmaxX = min_max_x(lines);
        return minmaxX.second - minmaxX.first;
    }

    double rangeY(Lines2D& lines){
        auto minmaxX = min_max_x(lines);
        return minmaxX.second - minmaxX.first;
    }

    std::pair<double, double> scaleLines(Lines2D& lines, const int& size){
        std::pair<double, double> minmax_x = min_max_x(lines);
        std::pair<double, double> minmax_y = min_max_y(lines);

        double biggestxy;
        if (minmax_x.second - minmax_x.first < minmax_y.second - minmax_y.first){
            biggestxy = minmax_y.second - minmax_y.first;
        }else{
            biggestxy = minmax_x.second - minmax_x.first;
        }

        double imagex = size*(minmax_x.second - minmax_x.first)/biggestxy;
        double imagey = size*(minmax_y.second - minmax_y.first)/biggestxy;

        double scalingfactor = .95*imagex/(minmax_x.second - minmax_x.first);

        double dc_x = scalingfactor *(minmax_x.first + minmax_x.second) / 2;
        double dc_y = scalingfactor *(minmax_y.first + minmax_y.second) / 2;

        double dx = imagex/2 - dc_x;
        double dy = imagey/2 - dc_y;

        for (auto& line : lines){
            line.p1.x *= scalingfactor; line.p1.x += dx; std::round(line.p1.x);
            line.p1.y *= scalingfactor; line.p1.y += dy; std::round(line.p1.y);
            line.p2.x *= scalingfactor; line.p2.x += dx; std::round(line.p2.x);
            line.p2.y *= scalingfactor; line.p2.y += dy; std::round(line.p2.y);
        }
        return std::pair<double, double>(imagex, imagey);
    }

    std::shared_ptr<img::EasyImage> draw2DLines(Lines2D &lines, const int size) {
        std::shared_ptr<img::EasyImage> image = std::make_shared<img::EasyImage>();
        scaleLines(lines, size);
        for( auto& i: lines){
            image->draw_line(i.p1.x, i.p1.y, i.p2.x, i.p2.y, i.color);
        }
        return image;
    }

} // namespace lines2D

