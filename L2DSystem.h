//
// Created by Phili on 2/26/2020.
//

#ifndef L2DSYSTEM_H_
#define L2DSYSTEM_H_


#include <fstream>
#include <iostream>
#include <memory>
#include<string>

#include "ini_configuration.h"
#include "l_parser.h"



namespace lsystem2D_gen {

    class Lsystem2D_gen: LParser::LSystem2D{
//        std::fstream L2D_file; //todo ifstream
//        LParser::LSystem2D lsys2d;
        public:
        explicit Lsystem2D_gen(std::istream& file) : LSystem2D(file) { file >> *this;}

        void GenLSys2D(const ini::Configuration& configuration);

        std::string getString();


    };

} // namespace lsystem2D_gen;


#endif //L2DSYSTEM_H_
