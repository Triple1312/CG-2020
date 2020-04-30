#include <fstream>
#include <iostream>
#include <string>
#include <chrono>

#include "easy_image.h"
#include "ini_configuration.h"

#include "l_parser.h"
#include "l_system.h"
#include "Figures3D.h"
#include "ZBuffer.h"



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
        //lines.resize(1023, lines_2d::Line2D());

        ini::DoubleTuple background = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        img::Color bg(uint8_t (background[0]*255), uint8_t (background[1]*255), uint8_t (background[2]*255));
//        bg.red = background[0]*256; bg.green = background[1]*256; bg.blue = background[2]*256;
        std::pair<int, int > size = scaleLines(lines, configuration["General"]["size"].as_int_or_die()); //todo voor hier
//        bg.red = uint8_t (background[0] *255);
        img::EasyImage image(size.first, size.second, bg);
        std::cout << lines.size() << std::endl;
        for (auto& i : lines){
            image.draw_line((unsigned int) tools::d2i(i.p1.x), (unsigned int) tools::d2i(i.p1.y),(unsigned int)tools::d2i(i.p2.x), (unsigned int)tools::d2i(i.p2.y), i.color);
        }
        return image;
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////

    else if(configuration["General"]["type"].as_string_or_die() == "Wireframe") {

        ini::DoubleTuple background = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
        img::Color bg(uint8_t (background[0]*255), uint8_t (background[1]*255), uint8_t (background[2]*255));
        lines_2d::Lines2D lines;

        figures_3d::Figures3D figs = figures_3d::Wireframe(configuration);

        for (auto &l : figs) {
            for (auto n : l.project_figure()){
                lines.push_back( n );
            }
        }

        std::pair<int, int > size = l_system2d::scaleLines(lines, configuration["General"]["size"].as_int_or_die());
        img::EasyImage image(size.first, size.second, bg);

        for (const auto& m : lines) {
            image.draw_line( tools::d2i(m.p1.x), tools::d2i(m.p1.y), tools::d2i(m.p2.x), tools::d2i(m.p2.y), m.color );
            image(3, 5) = bg;

        }
        return image;
    }

    else if (configuration["General"]["type"].as_string_or_die() == "ZBufferedWireframe"){
      ini::DoubleTuple background = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
      img::Color bg(uint8_t (background[0]*255), uint8_t (background[1]*255), uint8_t (background[2]*255));
      figures_3d::Figures3D figs = figures_3d::Wireframe(configuration);
      lines_2d::Lines2D lines;

      for (auto &l : figs) {
        for (auto n : l.project_figure()){
          lines.push_back( n );
        }
      }

      std::pair<int, int > size = l_system2d::scaleLines(lines, configuration["General"]["size"].as_int_or_die());

      img::EasyImage image(size.first, size.second, bg);
      zbuffer::ZBuffer buffer(size.first, size.second);
      return buffer.draw_lines(lines, size.first, size.second);

    }

    else if(configuration["General"]["type"].as_string_or_die() == "LightedZBuffering" ||
            configuration["General"]["type"].as_string_or_die() == "ZBuffering") { // todo was hier bezig
      ini::DoubleTuple background = configuration["General"]["backgroundcolor"].as_double_tuple_or_die();
      img::Color bg(uint8_t (background[0]*255), uint8_t (background[1]*255), uint8_t (background[2]*255));
      figures_3d::Figures3D figs = figures_3d::Wireframe(configuration);
      Lights3D lights = Lights3D(configuration);

      if (configuration["General"]["clipping"].as_bool_or_default(false)) {
        double near = configuration["General"]["dNear"].as_double_or_die();
        ini::DoubleTuple eye = configuration["General"]["eye"].as_double_tuple_or_die();
        Vector3D eyev = Vector3D::point(eye[0], eye[1], eye[2]);
        double r = std::pow(std::pow(eye[0], 2.0) + std::pow(eye[1], 2.0) + std::pow(eye[2], 2.0), 0.5);
        for (auto &f : figs){
          f.to_eye_clip(std::atan2(eye[1], eye[0]), std::acos(eye[2] / r),  eyev);
          f = figures_3d::clipping(f, configuration);

          for (auto &p : f.points){
            p.x = p.x * near / -p.z;
            p.y = p.y * near / -p.z;
          }
        }
      }

      lines_2d::Lines2D lines;
      for (auto &l : figs) {
        for (auto n : l.project_figure()){
          lines.push_back( n );
        }
      }
      if (lines.empty()) { std::cout << "geen figuren te tekenen" << std::endl;}
      std::tuple<double, double, double, double, double> d_dx_dy_image = lines_2d::d_dx_dy(lines, configuration["General"]["size"].as_int_or_die());
      double d, dx, dy;int imagex; int imagey;
      std::tie(d, dx, dy, imagex, imagey) = d_dx_dy_image;
      img::EasyImage image(imagex, imagey);
      zbuffer::ZBuffer buffer(imagex, imagey);

      for(auto &l : figs ){
        std::vector<figures_3d::Face> newfaces;
        for(auto k = 0; k < l.faces.size(); k++) {
          if (l.faces[k].size() > 3) {
            for (auto j : figures_3d::triangulate(l.faces[k])) {
              newfaces.push_back(j);
            }
          } else {
            newfaces.emplace_back(l.faces[k]);
          }
        }
        l.faces = newfaces;

        /// if no lights were given and a figure has a color than we will place it as the ambient reflection color and
        /// a white ambient light will be made here.
        if (lights.empty()){lights.emplace_back(Light());}

        for (int f = 0; f < l.faces.size(); ++f) {
          figures_3d::draw_triangle_light(l.points[l.faces[f][0]],
                                          l.points[l.faces[f][1]],
                                          l.points[l.faces[f][2]],
                                          dx,
                                          dy,
                                          d,
                                          image,
                                          buffer,
                                          l.ambientReflection,
                                          l.diffuseReflection,
                                          l.specularReflection,
                                          l.reflectionCoefficient,
                                          lights
                                           );
        }
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
