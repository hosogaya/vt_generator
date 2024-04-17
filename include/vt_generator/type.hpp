#pragma once

#include <vector>
#include <eigen3/Eigen/Core>
#include <ifopt/ipopt_solver.h>
#include <cppad/cppad.hpp>

namespace vt_generator
{
using Scalar = double;
using AdVector = std::vector<CppAD::AD<Scalar>>;
using AdScalar = CppAD::AD<Scalar>;
using Vector = Eigen::VectorXd;
using Vector2 = Eigen::Vector2d;
using Jacobian = ifopt::Component::Jacobian;
using VecBound = ifopt::Component::VecBound;
using Bound = ifopt::Bounds;

class Index
{
public:
    Index() {}
    ~Index() {}

    int acc(const int ind) const {return size_*ind + acc_;}
    int steer(const int ind) const {return size_*ind + steer_;}
    int n(const int ind) const {return size_*ind + n_;}
    int xi(const int ind) const {return size_*ind + xi_;}
    int vx(const int ind) const {return size_*ind + vx_;}
    int vy(const int ind) const {return size_*ind + vy_;}
    int w(const int ind) const {return size_*ind + w_;}
    int size() const {return size_;}
private:
    const int acc_ = 0;
    const int steer_ = 1;
    const int n_ = 2;
    const int xi_ = 3;
    const int vx_ = 4;
    const int vy_ = 5;
    const int w_ = 6;
    const int size_ = 7;

};

Index var_index;


}