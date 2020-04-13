//
// Created by Phili on 13/04/2020.
//

#ifndef ENGINE__ZBUFFER_H_
#define ENGINE__ZBUFFER_H_


#include <vector>
#include <limits>
#include <iostream>
#include <assert.h>
#include <utility>
#include <cmath>

#include "easy_image.h"
#include "vector3d.h"
#include "lines.h"

namespace zbuffer{
 class ZBuffer {

   std::vector<std::vector<double>> buffer;

  public:

   ZBuffer() = default;

   ZBuffer(int width, int height);

   img::EasyImage draw_lines(lines_2d::Lines2D &lines, const int&, const int&);

   double& operator()(int width, int height);

  };

}


#endif //ENGINE__ZBUFFER_H_
