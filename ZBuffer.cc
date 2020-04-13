//
// Created by Phili on 13/04/2020.
//

#include "ZBuffer.h"

zbuffer::ZBuffer::ZBuffer(int width, int height) {
  std::vector<std::vector<double>> temp(width, std::vector<double>(height, std::numeric_limits<double>::infinity()));
  this->buffer = temp;
}

img::EasyImage zbuffer::ZBuffer::draw_lines( lines_2d::Lines2D &lines, const int &width, const int &height ) {
  img::EasyImage image(width, height);
  int linecount = 0;

  for (auto l : lines){
    linecount +=1;
    if (linecount% 1000 == 0) {
      std::cout << linecount << std::endl << std::flush;
    }
    auto color = l.color;

    int x0 = std::round(l.p1.x);
    int y0 = std::round(l.p1.y);
    int z0 = std::round(l.z0);

    int x1 = std::round(l.p2.x);
    int y1 = std::round(l.p2.y);
    int z1 = std::round(l.z1);

    //de afstand tussen de 2 punten
    double r = std::pow(std::pow(x0 - x1, 2) + std::pow(y0 - y1, 2) + std::pow(z0 - z1, 2), 0.5);

    if (x0 > x1)
    {
      //flip points if x1>x0: we want x0 to have the lowest value
      std::swap(x0, x1);
      std::swap(y0, y1);
      std::swap(z0, z1);
    }

    double z_r = (z0 - z1) / r;

    double m = ((double) y1 - (double) y0) / ((double) x1 - (double) x0);
    if (-1.0 <= m && m <= 1.0)
    {
      for (unsigned int i = 0; i <= (x1 - x0); i++)
      {
        if (1/ (z0 + z_r *i) <buffer[x0 + i][ round(y0 + m * i)]){
          image((unsigned int) x0 + i, (unsigned int) round(y0 + m * i)) = color;
          buffer[x0 + i][std::round(y0 + m * i)] = 1/ (z0 + z_r *i);
        }
      }
    }
    else if (m > 1.0)
    {
      for (unsigned int i = 0; i <= (y1 - y0); i++)
      {
        if (1/ (z0 + z_r *i) < buffer[round(x0 + (i / m))][y0 + i]){
          image(round(x0 + (i / m)), y0 + i) = color;
          buffer[round(x0 + (i / m))][y0 + i] = 1/ (z0 + z_r *i);
        }
      }
    }
    else if (m < -1.0)
    {
      for (unsigned int i = 0; i <= (y0 - y1); i++)
      {
        if (1/ (z0 + z_r *i) < buffer[round(x0 - (i / m))][y0 - i]){
          image(round(x0 - (i / m)), y0 - i) = color;
          buffer[round(x0 - (i / m))][y0 - i] = 1/ (z0 + z_r *i);
        }
      }
    }
  }
  return image;
}

double &zbuffer::ZBuffer::operator()(int width, int height) {
  return this->buffer[width][height];
}


