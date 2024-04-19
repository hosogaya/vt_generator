#pragma once

#include <vt_generator/type.hpp>
#include <vt_generator/path_smoother/decoder.hpp>

namespace vt_generator
{
namespace path_smoother
{
class Cost: public ifopt::CostTerm
{
public:
    Cost(const int horizon)
    : ifopt::CostTerm("cost"), horizon_(horizon)
    {
        
    }
    ~Cost() {}

    Scalar GetCost() const override
    {
        auto x = GetVariables()->GetComponent("variables")->GetValues();

        std::vector<Scalar> ax(x.size());
        for (int i=0; i<x.size(); ++i) ax[i] = x(i);

        Scalar value = 0.0;
        for (int i=0; i<GetRows(); ++i)
        {
            auto fx = decorder_.decodeX(ax[(i-width_+horizon_)%horizon_], (i-width_+horizon_)%horizon_);
            auto fy = decorder_.decodeY(ax[(i-width_+horizon_)%horizon_], (i-width_+horizon_)%horizon_);
            auto cx = decorder_.decodeX(ax[i], i);
            auto cy = decorder_.decodeY(ax[i], i);
            auto lx = decorder_.decodeX(ax[(i+width_)%horizon_], (i+width_)%horizon_);
            auto ly = decorder_.decodeY(ax[(i+width_)%horizon_], (i+width_)%horizon_);

            auto area4 = std::pow(2*((cx - fx)*(ly - cy) - (cy - fy)*(lx - cx)), 2.0);
            auto denominater = std::pow(fx - cx, 2.0) + std::pow(fy - cy, 2.0)
                             + std::pow(lx - cx, 2.0) + std::pow(ly - cy, 2.0)
                             + std::pow(fx - lx, 2.0) + std::pow(fy - ly, 2.0);

            value += area4/denominater;
        }
        
        return value;
    }

    void FillJacobianBlock(std::string var_name, Jacobian& jac) const override
    {
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Scalar> ax(x.size());
        for (int i=0; i<x.size(); ++i) ax[i] = x(i);
        auto ad_jac = getADfun().Jacobian(ax);
        for (int i=0; i<x.size(); ++i)
        {
            jac.coeffRef(0, i) = ad_jac.at(i);
        }
    }

private:
    CppAD::ADFun<Scalar> getADfun() const 
    {
        AdVector ax(horizon_);
        AdVector ay(1);
        CppAD::Independent(ax);

        for (int i=0; i<horizon_; ++i)
        {
            auto fx = decorder_.decodeX(ax[(i-width_+horizon_)%horizon_], (i-width_+horizon_)%horizon_);
            auto fy = decorder_.decodeY(ax[(i-width_+horizon_)%horizon_], (i-width_+horizon_)%horizon_);
            auto cx = decorder_.decodeX(ax[i], i);
            auto cy = decorder_.decodeY(ax[i], i);
            auto lx = decorder_.decodeX(ax[(i+width_)%horizon_], (i+width_)%horizon_);
            auto ly = decorder_.decodeY(ax[(i+width_)%horizon_], (i+width_)%horizon_);

            auto area4 = CppAD::pow(2*((cx - fx)*(ly - cy) - (cy - fy)*(lx - cx)), 2.0);
            auto denominater = CppAD::pow(fx - cx, 2.0) + CppAD::pow(fy - cy, 2.0)
                             + CppAD::pow(lx - cx, 2.0) + CppAD::pow(ly - cy, 2.0)
                             + CppAD::pow(fx - lx, 2.0) + CppAD::pow(fy - ly, 2.0);

            ay[0] += area4/denominater;
        }

        return CppAD::ADFun<Scalar>(ax, ay);
    }
    
    const int horizon_;
    const int width_ = 3;
};
}
}