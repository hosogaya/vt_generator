#pragma once

#include <vt_generator/type.hpp>

namespace vt_generator
{
namespace path_smoother
{
class Decorder
{
public:
    Decorder() {}
    ~Decorder() {}

    Vector2 decode(const Scalar& x, const int index)
    {
        Vector2 pos;
        pos.x() = 0.5*((1.0 - x)*bounds_.at(2*index  ).lower_ + (1.0 + x)*bounds_.at(2*index  ).upper_);
        pos.y() = 0.5*((1.0 - x)*bounds_.at(2*index+1).lower_ + (1.0 + x)*bounds_.at(2*index+1).upper_);

        return pos;
    }

    template <typename T>
    T decodeX(const T& x, const int index)
    {
        return 0.5*((1.0 - x)*bounds_.at(2*index  ).lower_ + (1.0 + x)*bounds_.at(2*index  ).upper_);
    }
    template <typename T>
    T decodeY(const T& x, const int index)
    {
        return 0.5*((1.0 - x)*bounds_.at(2*index+1).lower_ + (1.0 + x)*bounds_.at(2*index+1).upper_);
    }

    void setBounds(const VecBound& b) 
    {
        bounds_ = b;
    }
private:    
    VecBound bounds_;
};

Decorder decorder_;
}
}