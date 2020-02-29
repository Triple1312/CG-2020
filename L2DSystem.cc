//
// Created by Phili on 2/26/2020.
//

#include "L2DSystem.h"

namespace lsystem2D_gen{

    void Lsystem2D_gen::GenLSys2D(const ini::Configuration& configuration) {
        std::fstream file;
        file.open(configuration["2DLSystem"]["inputfile"].as_string_or_die(),std::ios::in);
    }

    std::string Lsystem2D_gen::getString() {
        std::string string1 = this->get_initiator();
        std::string string2; // todo betere naam
        for (unsigned int i =0; i < this->get_nr_iterations();i++){
            for (const char j:string1){
                if(j == '+' or '-'){
                    string2 += j;
                }
                else{
                    string2 += this->get_replacement(j);
                }
            }
            string1 = string2;
            string2.clear();
        }
        return string1;
    }

}; //namespace lsystem2D_gen;

