#pragma once

#include <vt_generator/type.hpp>

namespace vt_generator
{
class TimeCost: public ifopt::CostTerm
{
public:
    TimeCost(const int horizon,
            const std::vector<Scalar>& curvature,
            const std::vector<Scalar>& ds)
    : ifopt::CostTerm("cost"), horizon_(horizon), curvature_(curvature), ds_(ds)
    {
        
    }
    ~TimeCost() {}

    Scalar GetCost() const override
    {
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Scalar> ax(x.size());
        for (int i=0; i<x.size(); ++i) 
        {
            ax[i] = x(i);
        }
        std::vector<Scalar> fw = getADFun().Forward(0, ax);
        assert(fw.size() == 1);
        return fw[0];
    }

    void FillJacobianBlock(std::string var_name, Jacobian& jac) const override
    {
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Scalar> ax(x.size());
        for (int i=0; i<x.size(); ++i)
        {
            ax[i] = x(i);
        }
        std::vector<Scalar> ad_jac = getADFun().Jacobian(ax);

        for (int i=0; i<horizon_; ++i)
        {
            jac.coeffRef(0, denormalizer_.acc(i))   = ad_jac[denormalizer_.acc(i)];
            jac.coeffRef(0, denormalizer_.steer(i)) = ad_jac[denormalizer_.steer(i)];
            jac.coeffRef(0, denormalizer_.n(i))     = ad_jac[denormalizer_.n(i)];
            jac.coeffRef(0, denormalizer_.xi(i))    = ad_jac[denormalizer_.xi(i)];
            jac.coeffRef(0, denormalizer_.vx(i))    = ad_jac[denormalizer_.vx(i)];
            jac.coeffRef(0, denormalizer_.vy(i))    = ad_jac[denormalizer_.vy(i)];
        }
    }

private:
    CppAD::ADFun<Scalar> getADFun() const
    {
        AdVector ay(1);
        AdVector ax(denormalizer_.size()*horizon_);
        CppAD::Independent(ax);

        ay[0] = 0;
        for (int i=0; i<horizon_; ++i)
        {
            auto n  = denormalizer_.denormalizeN(ax, i);
            auto xi = denormalizer_.denormalizeXi(ax, i);
            auto vx = denormalizer_.denormalizeVx(ax, i);
            auto vy = denormalizer_.denormalizeVy(ax, i);
            auto& k  = curvature_[i];

            auto vx_inv = 1.0/(vx + CppAD::log(1 + CppAD::exp(-2*vx*alpha_))/alpha_);
            auto v = CppAD::pow(vx, 2) + CppAD::pow(vy, 2);
            auto beta = CppAD::atan(vy*vx_inv);
            auto ds = (v*CppAD::cos(xi + beta));
            auto dtds = (1.0 - n*k) / (ds + CppAD::log(1 + CppAD::exp(-2*ds*alpha_))/alpha_);

            // auto v = CppAD::pow(vx, 2) + CppAD::pow(vy, 2);
            // auto beta = CppAD::atan2(vy, vx);
            // auto dtds = (1.0 - n*k) / (v*CppAD::cos(xi + beta));
            ay[0] += dtds*ds_[i];
            ay[0] += CppAD::pow(ax.at(denormalizer_.acc((i+1)%horizon_))   - ax.at(denormalizer_.acc(i)),   2.0);
            ay[0] += CppAD::pow(ax.at(denormalizer_.steer((i+1)%horizon_)) - ax.at(denormalizer_.steer(i)), 2.0);
        }
        return CppAD::ADFun<Scalar>(ax, ay);
    }

    const int horizon_;
    std::vector<Scalar> ds_;
    std::vector<Scalar> curvature_;
    const Scalar alpha_ = 1e3;
};
}