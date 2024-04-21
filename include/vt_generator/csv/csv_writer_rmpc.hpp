#pragma once

#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <vt_generator/type.hpp>

namespace vt_generator
{
namespace csv
{
class WriterRMPC
{
public:
    WriterRMPC(const std::string& file_name);
    ~WriterRMPC() {}

    bool isOpen() const {return is_open_;}
    void writeResult(const std::vector<Scalar>& outer_w,  const std::vector<Scalar>& inner_w,
                     const std::vector<Scalar>& outer_x,  const std::vector<Scalar>& outer_y,
                     const std::vector<Scalar>& inner_x,  const std::vector<Scalar>& inner_y,
                     const std::vector<std::vector<Scalar>>& opt_state,
                     const std::vector<std::pair<Scalar, Scalar>>& opt_input,
                     const std::vector<Scalar>& curvature);
private:
    bool is_open_ = false;
    std::ofstream output_file_;
};

WriterRMPC::WriterRMPC(const std::string& file_name)
{
    output_file_.open(file_name);
    if (!output_file_) 
    {
        is_open_ = false;
        return;
    }
    is_open_ = true;
}

// opt_x, opt_y, outer_width, inner_width
// center_x, center_y, outer_x, outer_y inner_x, inner_y
// curvature, ref_v
void WriterRMPC::writeResult(const std::vector<Scalar>& outer_w,  const std::vector<Scalar>& inner_w,
                        const std::vector<Scalar>& outer_x,  const std::vector<Scalar>& outer_y,
                        const std::vector<Scalar>& inner_x,  const std::vector<Scalar>& inner_y,
                        const std::vector<std::vector<Scalar>>& opt_state,
                        const std::vector<std::pair<Scalar, Scalar>>& opt_input,
                        const std::vector<Scalar>& curvature)
{
    output_file_ << "outer_width" << ","
                 << "inner_width" << ","
                 << "outer_x"     << ","
                 << "outer_y"     << ","
                 << "inner_x"     << ","
                 << "inner_y"     << ","
                 << "acc"     << ","
                 << "steer"     << ","
                 << "normal"     << ","
                 << "xi"     << ","
                 << "vx"     << ","
                 << "vy"     << ","
                 << "w"     << ","
                 << "curvature"  << std::endl;

    for (int i=0; i<outer_w.size(); ++i)
    {
        output_file_ << outer_w[i]   << ","
                     << inner_w[i]   << ","
                     << outer_x[i]   << ","
                     << outer_y[i]   << ","
                     << inner_x[i]   << ","
                     << inner_y[i]   << ","
                     << opt_input[i].first << ","
                     << opt_input[i].second << ","
                     << opt_state[i][0] << ","
                     << opt_state[i][1] << ","
                     << opt_state[i][2] << ","
                     << opt_state[i][3] << ","
                     << opt_state[i][4] << ","
                     << curvature[i] << std::endl;
    }
}

}
}