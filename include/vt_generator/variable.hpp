#pragma once

#include <vt_generator/type.hpp>

namespace vt_generator
{
class Variable: public ifopt::VariableSet
{
public: 
    Variable(const int horizon, const VecBound& b) :  
        horizon_(horizon), bounds_(b), VariableSet(var_index.size()*horizon, "variables")
    {
		    // acc, steer
        // n, xi, vx, vy, omega
        x_.resize(var_index.size()*horizon);
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
    const VecBound bounds_;
};

}