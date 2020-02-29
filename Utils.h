//
// Created by Phili on 2/19/2020.
//

#ifndef ENGINE_UTILS_H
#define ENGINE_UTILS_H


#include <iostream>
#include "easy_image.h"


struct Color{
    double red;
    double green;
    double blue;
};

struct Point2D{
    double x;
    double y;
};

struct Line2D{
    Point2D p1;
    Point2D p2;
    Color color;
};

using Lines2D = std::vector<Line2D>;

double rangeX(Lines2D lines){
    double xmin = lines[0].p1.x;
    double xmax = lines[0].p1.x;
    for (auto xp : lines){  // kan verbeterd worden
        if (xp.p1.x < xmin){ xmin = xp.p1.x;}
        if (xp.p1.x > xmax){ xmax = xp.p1.x;}
        if (xp.p2.x < xmin){ xmin = xp.p2.x;}
        if (xp.p2.x > xmax){ xmax = xp.p2.x;}
    }
    return xmax - xmin;
}

double rangeY(Lines2D lines){
    double ymin = lines[0].p1.y;
    double ymax = lines[0].p1.y;
    for (auto yp : lines){
        if (yp.p1.y < ymin){ ymin = yp.p1.y;}
        if (yp.p1.y > ymax){ ymax = yp.p1.y;}
        if (yp.p2.y < ymin){ ymin = yp.p2.y;}
        if (yp.p2.y > ymax){ ymax = yp.p2.y;}
    }
    return ymax - ymin;
}

void scaleLines(Lines2D lines, int size){
    double rangex = rangeX(lines);
    double rangey = rangeY(lines);
    double biggestxy;
    if (rangex < rangey){
        biggestxy = rangey;
    }else{biggestxy = rangex;
    }

    double imagex = size*rangex/biggestxy;
    double imagey = size*rangey/biggestxy;

    double scalingfactor = .95*imagex/rangex;

    for (auto line : lines){
        line.p1.x *= scalingfactor;
        line.p1.y *= scalingfactor;
        line.p2.x *= scalingfactor;
        line.p2.y *= scalingfactor;
    }
}

img::EasyImage draw2DLines(const Lines2D& lines,const int size){
    double rangex = rangeX(lines);
    double rangey = rangeY(lines);
    double biggestxy;
    if (rangex < rangey){
        biggestxy = rangey;
    }else{biggestxy = rangex;
    }

    double imagex = size*rangex/biggestxy;
    double imagey = size*rangey/biggestxy;

    double scalingfactor = .95*imagex/rangex;

    for (auto line : lines){
        line.p1.x *= scalingfactor;
        line.p1.y *= scalingfactor;
        line.p2.x *= scalingfactor;
        line.p2.y *= scalingfactor;
    }
     double DCx =




}



#endif //ENGINE_UTILS_H
