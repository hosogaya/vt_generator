#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

namespace vt_generator
{
namespace csv
{
class Reader
{
public:
    Reader(const std::string& file_name, bool has_header = true);
    ~Reader() {}

    const std::vector<double>& xm() const {return xm_;}
    const std::vector<double>& ym() const {return ym_;}
    const std::vector<double>& outer_x() const {return outer_x_;}
    const std::vector<double>& outer_y() const {return outer_y_;}
    const std::vector<double>& inner_x() const {return inner_x_;}
    const std::vector<double>& inner_y() const {return inner_y_;}
    const std::vector<double>& curvature() const {return curvature_;}
    const std::vector<double>& outer_w() const {return outer_width_;}
    const std::vector<double>& inner_w() const {return inner_width_;}
    const std::vector<double>& center_x() const {return center_x_;}
    const std::vector<double>& center_y() const {return center_y_;}
    const std::vector<double>& ref_v() const {return ref_v_;}
    int size() const {return size_;}
    bool isOpen() const {return is_open_;}

private:
    std::vector<double> xm_;
    std::vector<double> ym_;
    std::vector<double> inner_width_;
    std::vector<double> outer_width_;
    std::vector<double> center_x_;
    std::vector<double> center_y_;
    std::vector<double> outer_x_;
    std::vector<double> outer_y_;
    std::vector<double> inner_x_;
    std::vector<double> inner_y_;
    std::vector<double> curvature_;
    std::vector<double> ref_v_;
    int size_;
    bool is_open_ = false;
};

Reader::Reader(const std::string& file_name, bool has_header)
{
    std::ifstream ifs(file_name);
    if (!ifs.is_open())
    {
        std::cout << "Failed to open " << file_name << std::endl;
        return;
    }
    std::cout << "Open " << file_name << std::endl;
    is_open_ = true;
    std::string str_buf;
    std::string str_conma_buf;
    size_ = 0;
    if (has_header) std::getline(ifs, str_buf); // skip header
    while (std::getline(ifs, str_buf))
    {
        ++size_;
        std::istringstream i_stream(str_buf);
        // 
        std::getline(i_stream, str_conma_buf, ',');
        xm_.push_back(std::stod(str_conma_buf));
        
        std::getline(i_stream, str_conma_buf, ',');
        ym_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf, ',');
        outer_width_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf, ',');
        inner_width_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf, ',');
        center_x_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf, ',');
        center_y_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf, ',');
        outer_x_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf, ',');
        outer_y_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf, ',');
        inner_x_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf, ',');
        inner_y_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf, ',');
        curvature_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf, ',');
        ref_v_.push_back(std::stod(str_conma_buf));
    }
}

}
}