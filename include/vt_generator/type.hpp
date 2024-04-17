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

class Denormalizer
{
public:
    // Index() {}
    // ~Index() {}

    Denormalizer() {}
    ~Denormalizer() {}

    int acc(const int ind) const {return size_*ind + acc_;}
    int steer(const int ind) const {return size_*ind + steer_;}
    int n(const int ind) const {return size_*ind + n_;}
    int xi(const int ind) const {return size_*ind + xi_;}
    int vx(const int ind) const {return size_*ind + vx_;}
    int vy(const int ind) const {return size_*ind + vy_;}
    int w(const int ind) const {return size_*ind + w_;}
    int size() const {return size_;}

    void setVariableBound(const VecBound& b) {bounds_ = b;}
    template <typename T> 
    T denormalizeAcc(const std::vector<T>& x, const int ind)
    {
        return denormalize(x.at(acc(ind)), bounds_.at(acc(ind)));
    }

    template <typename T> 
    T denormalizeSteer(const std::vector<T>& x, const int ind)
    {
        return denormalize(x.at(steer(ind)), bounds_.at(steer(ind)));
    }

    template <typename T> 
    T denormalizeN(const std::vector<T>& x, const int ind)
    {
        return denormalize(x.at(n(ind)), bounds_.at(n(ind)));
    }

    template <typename T> 
    T denormalizeXi(const std::vector<T>& x, const int ind)
    {
        return denormalize(x.at(xi(ind)), bounds_.at(xi(ind)));
    }

    template <typename T> 
    T denormalizeVx(const std::vector<T>& x, const int ind)
    {
        return denormalize(x.at(vx(ind)), bounds_.at(vx(ind)));
    }

    template <typename T> 
    T denormalizeVy(const std::vector<T>& x, const int ind)
    {
        return denormalize(x.at(vy(ind)), bounds_.at(vy(ind)));
    }

    template <typename T> 
    T denormalizeW(const std::vector<T>& x, const int ind)
    {
        return denormalize(x.at(w(ind)), bounds_.at(w(ind)));
    }

    template <typename T>
    T denormalize(const T& x, const Bound& b)
    {
        return 0.5*((1.0 - x)*b.lower_ + (1.0 + x)*b.upper_);
    }
private:
    const int acc_ = 0;
    const int steer_ = 1;
    const int n_ = 2;
    const int xi_ = 3;
    const int vx_ = 4;
    const int vy_ = 5;
    const int w_ = 6;
    const int size_ = 7;

    VecBound bounds_;
};

Denormalizer denormalizer_;


const int constraints_size_ = 7;

}