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
        std::cout << "get value of cost" << std::endl;
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
        std::cout << "fill jacobian of cost" << std::endl;
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Scalar> ax(x.size());
        for (int i=0; i<x.size(); ++i)
        {
            ax[i] = x(i);
        }
        std::vector<Scalar> ad_jac = getADFun().Jacobian(ax);

        for (int i=0; i<horizon_; ++i)
        {
            jac.coeffRef(0, var_index.n(i))  = ad_jac[var_index.n(i)];
            jac.coeffRef(0, var_index.xi(i)) = ad_jac[var_index.xi(i)];
            jac.coeffRef(0, var_index.vx(i)) = ad_jac[var_index.vx(i)];
            jac.coeffRef(0, var_index.vy(i)) = ad_jac[var_index.vy(i)];
        }
    }

private:
    CppAD::ADFun<Scalar> getADFun() const
    {
        AdVector ay(1);
        AdVector ax(var_index.size()*horizon_);
        CppAD::Independent(ax);

        ay[0] = 0;
        for (int i=0; i<horizon_; ++i)
        {
            auto& n = ax[var_index.n(i)];
            auto& xi = ax[var_index.xi(i)];
            auto& vx = ax[var_index.vx(i)];
            auto& vy = ax[var_index.vy(i)];
            auto& k = curvature_[i];

            auto vx_inv = 1.0/(vx + CppAD::log(1 + CppAD::exp(-2*vx*alpha_))/alpha_);
            auto v = CppAD::pow(vx, 2) + CppAD::pow(vy, 2);
            auto beta = CppAD::atan(vy*vx_inv);
            auto ds = (v*CppAD::cos(xi + beta));
            auto dtds = (1.0 - n*k) / (ds + CppAD::log(1 + CppAD::exp(-2*ds*alpha_))/alpha_);

            // auto v = CppAD::pow(vx, 2) + CppAD::pow(vy, 2);
            // auto beta = CppAD::atan2(vy, vx);
            // auto dtds = (1.0 - n*k) / (v*CppAD::cos(xi + beta));
            ay[0] += dtds*ds_[i];
        }
        return CppAD::ADFun<Scalar>(ax, ay);
    }

    const int horizon_;
    std::vector<Scalar> ds_;
    std::vector<Scalar> curvature_;
    const Scalar alpha_ = 1e3;
};
}