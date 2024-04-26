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
class Writer
{
public:
    Writer(const std::string& file_name);
    ~Writer() {}

    bool isOpen() const {return is_open_;}
    void writeResult(const Vector& x_opt, const int horizon);
    void writeResult(const std::vector<Scalar>& opt_x,    const std::vector<Scalar>& opt_y,
                     const std::vector<Scalar>& outer_w,  const std::vector<Scalar>& inner_w,
                     const std::vector<Scalar>& center_x, const std::vector<Scalar>& center_y,
                     const std::vector<Scalar>& outer_x,  const std::vector<Scalar>& outer_y,
                     const std::vector<Scalar>& inner_x,  const std::vector<Scalar>& inner_y,
                     const std::vector<Scalar>& curvature,const std::vector<Scalar>& ref_v,
                     const std::vector<Scalar>& true_outer_x, const std::vector<Scalar>& true_outer_y,
                     const std::vector<Scalar>& true_inner_x, const std::vector<Scalar>& true_inner_y,
                     const int horizon);
private:
    bool is_open_ = false;
    std::ofstream output_file_;
};

Writer::Writer(const std::string& file_name)
{
    output_file_.open(file_name);
    if (!output_file_) 
    {
        is_open_ = false;
        return;
    }
    is_open_ = true;
}

void Writer::writeResult(const Vector& x_opt, const int horizon)
{
    output_file_ << "acc"   << ","
                 << "steer" << ","
                 << "n"     << ","
                 << "xi"    << ","
                 << "vx"    << ","
                 << "vy"    << ","
                 << "w"     << std::endl;
    std::vector<Scalar> ax(x_opt.size());
    for (int i=0; i<x_opt.size(); ++i) ax[i] = x_opt(i);
    for (int i=0; i<horizon; ++i)
    {
        output_file_ << denormalizer_.denormalizeAcc(ax, i)  << ",";
        output_file_ << denormalizer_.denormalizeSteer(ax, i) << ",";
        output_file_ << denormalizer_.denormalizeN(ax, i)    << ",";
        output_file_ << denormalizer_.denormalizeXi(ax, i)   << ",";
        output_file_ << denormalizer_.denormalizeVx(ax, i)   << ",";
        output_file_ << denormalizer_.denormalizeVy(ax, i)   << ",";
        output_file_ << denormalizer_.denormalizeW(ax, i)    << std::endl;
    }
}

// opt_x, opt_y, outer_width, inner_width
// center_x, center_y, outer_x, outer_y inner_x, inner_y
// curvature, ref_v
void Writer::writeResult(const std::vector<Scalar>& opt_x,    const std::vector<Scalar>& opt_y,
                         const std::vector<Scalar>& outer_w,  const std::vector<Scalar>& inner_w,
                         const std::vector<Scalar>& center_x, const std::vector<Scalar>& center_y,
                         const std::vector<Scalar>& outer_x,  const std::vector<Scalar>& outer_y,
                         const std::vector<Scalar>& inner_x,  const std::vector<Scalar>& inner_y,
                         const std::vector<Scalar>& curvature,const std::vector<Scalar>& ref_v, 
                         const std::vector<Scalar>& true_outer_x, const std::vector<Scalar>& true_outer_y,
                         const std::vector<Scalar>& true_inner_x, const std::vector<Scalar>& true_inner_y,
                         const int horizon)
{
    output_file_ << "opt_x"       << ","
                 << "opt_y"       << ","
                 << "outer_width" << ","
                 << "inner_width" << ","
                 << "center_x"    << ","
                 << "center_y"    << ","
                 << "outer_x"     << ","
                 << "outer_y"     << ","
                 << "inner_x"     << ","
                 << "inner_y"     << ","
                 << "curvature"   << ","
                 << "ref_v"       << ","
                 << "true_outer_x"       << ","
                 << "true_outer_y"       << ","
                 << "true_inner_x"       << ","
                 << "true_inner_y"       << std::endl;

    for (int i=0; i<horizon; ++i)
    {
        output_file_ << opt_x[i]     << ","
                     << opt_y[i]     << ","
                     << outer_w[i]   << ","
                     << inner_w[i]   << ","
                     << center_x[i]  << ","
                     << center_y[i]  << ","
                     << outer_x[i]   << ","
                     << outer_y[i]   << ","
                     << inner_x[i]   << ","
                     << inner_y[i]   << ","
                     << curvature[i] << ","
                     << ref_v[i]     << ","
                     << true_outer_x[i]     << ","
                     << true_outer_y[i]     << ","
                     << true_inner_x[i]     << ","
                     << true_inner_y[i]     << std::endl;
    }
}

}
}