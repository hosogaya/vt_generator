#pragma once

#include <fstream>
#include <sstream>
#include <vector>

namespace vt_generator
{
namespace csv
{
class Reader
{
public:
    Reader(const std::string& file_name);
    ~Reader() {}

    const std::vector<double>& xm() const {return xm_;}
    const std::vector<double>& ym() const {return ym_;}
    const std::vector<double>& tr_right() const {return tr_right_;}
    const std::vector<double>& tr_left() const {return tr_left_;}
    int size() const {return size_;}

private:
    std::vector<double> xm_;
    std::vector<double> ym_;
    std::vector<double> tr_left_;
    std::vector<double> tr_right_;
    int size_;
};

Reader::Reader(const std::string& file_name)
{
    std::ifstream ifs(file_name);
    std::string str_buf;
    std::string str_conma_buf;
    size_ = 0;
    while (std::getline(ifs, str_buf))
    {
        ++size_;
        std::istringstream i_stream(str_buf);
        // 
        std::getline(i_stream, str_conma_buf);
        xm_.push_back(std::stod(str_conma_buf));
        
        std::getline(i_stream, str_conma_buf);
        ym_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf);
        tr_right_.push_back(std::stod(str_conma_buf));

        std::getline(i_stream, str_conma_buf);
        tr_left_.push_back(std::stod(str_conma_buf));
    }
}

}
}