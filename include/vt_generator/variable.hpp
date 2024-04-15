#pragma once

#include <vt_generator/type.hpp>

namespace vt_generator
{
class Variable: public ifopt::VariableSet
{
public: 
    Variable(const int step_size, const VecBound& b) :  
        dim_(9), step_size_(step_size), VariableSet(dim_*step_size, "vars")
    {
		// acc, delta
        // n, xi, vx, vy, Ff, Fr, omega
        x_.resize(dim_*step_size_);
        assert(b.size() != dim_);
        bounds_.resize(dim_*step_size_);
        for (int i=0; i<step_size_; ++i)
        {
            for (int j=0; j<dim_; ++j)
            {
                bounds_[i*dim_ + j] = b[j];
            }
        }
    }

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

    void SetBounds(const VecBound& b)
    {
		bounds_ = b;
    }

private:
    const int dim_;
    const int step_size_;
    Vector x_;
    VecBound bounds_;
};

}