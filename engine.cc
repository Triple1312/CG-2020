#include <fstream>
#include <iostream>
#include <string>
#include <chrono>

#include "easy_image.h"
#include "ini_configuration.h"

#include "l_parser.h"
#include "l_system.h"
#include "Figures3D.h"



///based on https://www.youtube.com/watch?v=YG4jexlSAjc&list=PLlrATfBNZ98dudnM48yfGUldqGD0S4FFb&index=74
class Timer{

public:
    Timer(){
        m_startpoint = std::chrono::high_resolution_clock::now();
    }

     ~Timer(){
        Stop();
    }

    void Stop(){
        auto endpoint = std::chrono::high_resolution_clock::now();

        auto start = std::chrono::time_point_cast<std::chrono::milliseconds>(m_startpoint).time_since_epoch().count();
        auto end = std::chrono::time_point_cast<std::chrono::milliseconds>(endpoint).time_since_epoch().count();
        std::cout << end - start;
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> m_startpoint;

};

img::EasyImage generate_image(const ini::Configuration &configuration) {
    Timer timer;

    if(configuration["General"]["type"].as_string_or_die() == "2DLSystem"){
        using namespace::l_system2d;

        std::vector<int> blub;

        LParser::LSystem2D l_system;
        std::ifstream file(configuration["2DLSystem"]["inputfile"].as_string_or_die());
        file >> l_system;
        ini::DoubleTuple lineColor = configuration["2DLSystem"]["color"].as_double_tuple_or_die();
        img::Color line_c(lineColor[0]*255, lineColor[1]*255, lineColor[2]*255);
        Lines2D lines = calcLSystem(l_system, line_c);
        lines.resize(1024, lines_2d::Line2D());

        ini::DoubleTuple background = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        img::Color bg(uint8_t (background[0]*255), uint8_t (background[1]*255), uint8_t (background[2]*255));
//        bg.red = background[0]*256; bg.green = background[1]*256; bg.blue = background[2]*256;
        std::pair<int, int > size = scaleLines(lines, configuration["General"]["size"].as_int_or_die());
//        bg.red = uint8_t (background[0] *255);
        img::EasyImage image(size.first, size.second, bg);
        for (auto& i : lines){
            image.draw_line(i.p1.x, i.p1.y, i.p2.x, i.p2.y, i.color);
        }
        return image;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    else if(configuration["General"]["type"].as_string_or_die() == "Wireframe") {

        figures_3d::Figures3D figs;
        int blub = configuration["General"]["nrFigures"].as_int_or_die();

        ini::DoubleTuple eye = configuration["General"]["eye"].as_double_tuple_or_die();
        for (auto i = 0; i < configuration["General"]["nrFigures"].as_int_or_die(); i++){

            figures_3d::Figure figure;  // todo kan veel beter door direct to initializen in de vector
            ini:: DoubleTuple fig_color = configuration["Figure" + std::to_string(i)]["color"].as_double_tuple_or_die();
            figure.color = img::Color(fig_color[0]*255, fig_color[1]*255, fig_color[2]*255);

            for (auto j=0; j < configuration["Figure" + std::to_string(i)]["nrPoints"].as_int_or_die() ; j++){
                ini::DoubleTuple point = configuration["Figure" + std::to_string(i)]["point" + std::to_string(j)].as_double_tuple_or_die();
                //Vector3D newPoint = Vector3D::point(point[0], point[1], point[2]);
                figure.points.push_back(Vector3D::point(point[0], point[1], point[2]));
            }

            for (auto k = 0; k < configuration["Figure" + std::to_string(i)]["nrLines"].as_int_or_die(); k++ ) {
                ini::IntTuple line = configuration["Figure" + std::to_string(i)]["line" + std::to_string(k)].as_int_tuple_or_die();
                //Vector3D newPoint = Vector3D::point(point[0], point[1], point[2]);
                std::vector<int> line_n; line_n.push_back(line[0]); line_n.push_back(line[1]);
                figure.faces.push_back(line_n);
            }
            double r = std::pow(std::pow(eye[0], 2.0) + std::pow(eye[1], 2.0) + std::pow(eye[2], 2.0), 0.5);
            figure.scale(configuration["Figure" + std::to_string(i)]["scale"].as_double_or_die());
            figure.rotateX(M_PI /180 * configuration["Figure" + std::to_string(i)]["rotateX"].as_int_or_die());
            figure.rotateY(M_PI /180 * configuration["Figure" + std::to_string(i)]["rotateY"].as_int_or_die());
            figure.rotateZ(M_PI /180 * configuration["Figure" + std::to_string(i)]["rotateZ"].as_int_or_die());
            figure.to_eye( std::atan2(eye[1],eye[0]), std::acos(eye[2]/ r), r);

            figs.push_back(figure);
        }

        ini::DoubleTuple background = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        img::Color bg(uint8_t (background[0]*255), uint8_t (background[1]*255), uint8_t (background[2]*255));
        lines_2d::Lines2D lines;

        for (auto &l : figs) {
            for (auto n : l.project_figure()){
                lines.push_back( n );
            }
        }

        std::pair<int, int > size = l_system2d::scaleLines(lines, configuration["General"]["size"].as_int_or_die());
        img::EasyImage image(size.first, size.second, bg);

        for (auto m : lines) {
            image.draw_line( std::round(m.p1.x), std::round(m.p1.y), std::round(m.p2.x), std::round(m.p2.y), m.color );
        }
        return image;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////
    img::EasyImage blub;
    return blub;
}

int main(int argc, char const* argv[])
{
        int retVal = 0;
        try
        {
                for(int i = 1; i < argc; ++i)
                {
                        ini::Configuration conf;
                        try
                        {
                                std::ifstream fin(argv[i]);
                                fin >> conf;
                                fin.close();
                        }
                        catch(ini::ParseException& ex)
                        {
                                std::cerr << "Error parsing file: " << argv[i] << ": " << ex.what() << std::endl;
                                retVal = 1;
                                continue;
                        }

                        img::EasyImage image = generate_image(conf);
                        if(image.get_height() > 0 && image.get_width() > 0)
                        {
                                std::string fileName(argv[i]);
                                std::string::size_type pos = fileName.rfind('.');
                                if(pos == std::string::npos)
                                {
                                        //filename does not contain a '.' --> append a '.bmp' suffix
                                        fileName += ".bmp";
                                }
                                else
                                {
                                        fileName = fileName.substr(0,pos) + ".bmp";
                                }
                                try
                                {
                                        std::ofstream f_out(fileName.c_str(),std::ios::trunc | std::ios::out | std::ios::binary);
                                        f_out << image;

                                }
                                catch(std::exception& ex)
                                {
                                        std::cerr << "Failed to write image to file: " << ex.what() << std::endl;
                                        retVal = 1;
                                }
                        }
                        else
                        {
                                std::cout << "Could not generate image for " << argv[i] << std::endl;
                        }
                }
        }
        catch(const std::bad_alloc &exception)
        {
    		//When you run out of memory this exception is thrown. When this happens the return value of the program MUST be '100'.
    		//Basically this return value tells our automated test scripts to run your engine on a pc with more memory.
    		//(Unless of course you are already consuming the maximum allowed amount of memory)
    		//If your engine does NOT adhere to this requirement you risk losing points because then our scripts will
		//mark the test as failed while in reality it just needed a bit more memory
                std::cerr << "Error: insufficient memory" << std::endl;
                retVal = 100;
        }
        return retVal;
}
