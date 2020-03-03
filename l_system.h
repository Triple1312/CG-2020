//
// Created by Phili on 2/27/2020.
//

#ifndef L_SYSTEM_H_
#define L_SYSTEM_H_

#include <iostream>
#include <cmath>
#include <memory>
#include <stack>
#include <utility>
#include <string>
#include <tuple>

#include "l_parser.h"
#include "lines.h"

/// a namespace for 2d LSystems
namespace l_system2d {
    using namespace lines_2d;

    /// Function a
    /// \param l_system
    /// \param linecolor
    /// \return all the lines of the LSystem
   std::shared_ptr<Lines2D> calcLSystem(const LParser::LSystem2D &l_system, const img::Color& linecolor);

   /// Function that calculates the endpoints of a line
   /// \return the end x and end y coordinate of a line
   std::pair<double,double> l_calc_line(const double&, const double&, const double&, const double&);

}  // namespace l_system2d

namespace l_system3d {

}  // namespace l_system3d

#endif  // L_SYSTEM_H_
