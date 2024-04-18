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

    output_file_ << "acc"   << ","
                 << "steer" << ","
                 << "n"     << ","
                 << "xi"    << ","
                 << "vx"    << ","
                 << "vy"    << ","
                 << "w"     << std::endl;
}

void Writer::writeResult(const Vector& x_opt, const int horizon)
{
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

}
}