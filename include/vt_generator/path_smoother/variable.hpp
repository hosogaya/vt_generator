#pragma once

#include <vt_generator/type.hpp>
#include <vt_generator/path_smoother/decoder.hpp>

namespace vt_generator
{
namespace path_smoother
{
class Variable: public ifopt::VariableSet
{
public: 
    Variable(const int horizon) :  
        horizon_(horizon), VariableSet(horizon, "variable")
    {
        // acc, steer
        // n, xi, vx, vy, omega
        x_.resize(GetRows());
        bounds_.resize(GetRows());
        for (auto& bound: bounds_)
        {
            bound.lower_ =-1.0;
            bound.upper_ = 1.0;
        }
    }

    ~Variable() {}

    void SetVariables(const VectorXd& x) override
    {
        assert(x_.size() == x.size());
        x_ = x;
    }

    VectorXd GetValues() const override
    {
	  	return x_;
    }

    VecBound GetBounds() const override
    {
		return bounds_;
    }

private:
    const int horizon_;
    Vector x_;
    VecBound bounds_;
};
}
}