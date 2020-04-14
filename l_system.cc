//
// Created by Phili on 2/27/2020.
//

#include "l_system.h"

namespace {

    std::string getString(const LParser::LSystem& l_sys) {
        std::string string1 = l_sys.get_initiator();
        std::string string2;  // todo betere naam
        for (unsigned int i =0; i < l_sys.get_nr_iterations(); i++) {
            for (const char j : string1) {
                if (j == '+' || j == '-' || j == '(' || j == ')') {
                    string2 += j;
                } else {
                    string2 += l_sys.get_replacement(j);
                }
            }
            string1 = string2;
            string2.clear();
        }
        return string1;
    }

}  // namespace

namespace l_system2d {

    Lines2D calcLSystem(
            const LParser::LSystem2D& l_system, const img::Color& linecolor) {
        std::string full_function = getString(l_system);
        Lines2D lines  ;
        std::stack<std::tuple<double, double, double>> stack;
        double curr_angle = l_system.get_starting_angle() * M_PI /180;
        double curr_x(0.0), curr_y(0.0);
        double angle = l_system.get_angle()* M_PI / 180;
        int test = 0;

        int counter = 0;
        for (char& i : full_function) {

            if (i == '+') {
                curr_angle += angle;

            } else if (i == '-') {
                curr_angle -= angle;

            } else if (i == '(') {
                stack.push(std::tuple<double, double, double>
                        (curr_x, curr_y, curr_angle));

            } else if (i == ')') {
                std::tie(curr_x , curr_y, curr_angle) = stack.top();
                stack.pop();

            } else {
              test += 1;
                if (l_system.draw(i)) {
                    counter += 1;
                    auto end = l_calc_line(curr_x, curr_y, 1.0, curr_angle);
                    lines.push_back(Line2D(Point2D(curr_x, curr_y),
                                    Point2D(end), linecolor));

                    curr_x =  end.first;
                    curr_y = end.second;
                }
            }
        }
        std::cout << test << std::endl;
        std::cout << counter << std::endl;
        return lines;
    }

    inline std::pair<double, double> l_calc_line(
            const double& start_x, const double& start_y,
            const double& length, const double& angle) {
        return std::pair<double, double>(
                start_x + length*(cos(angle)), start_y + length*(sin(angle)));
    }
}  // namespace l_system2d
