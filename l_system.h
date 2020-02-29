//
// Created by Phili on 2/27/2020.
//

#ifndef ENGINE_L_SYSTEM_H
#define ENGINE_L_SYSTEM_H

#include <iostream>
#include <cmath>
#include <memory>
#include <stack>

#include "l_parser.h"
#include "lines.h"

//namespace{
//
//} // namespace

namespace l_system2d {
    using namespace lines_2d;

   std::shared_ptr<Lines2D> calcLSystem(const LParser::LSystem2D &l_system, img::Color& linecolor);

   std::pair<double,double> l_calc_line(const double&, const double&, const double&, const double&);



} // namespace l_system2d

namespace l_system3d{

} // namespace l_system3d

#endif //ENGINE_L_SYSTEM_H
