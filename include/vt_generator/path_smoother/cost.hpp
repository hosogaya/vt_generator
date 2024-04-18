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
        // Scalar value = 0.0;
        // std::vector<Vector2> ax(GetRows());
        // for (int i=0; i<horizon_; ++i) ax[i] = decorder_.decode(x(i), i);
        // for (int i=0; i<horizon_; ++i)
        // {
        //    int former = (i-1+horizon_)%horizon_;
        //    int later  = (i+1)%horizon_;
        //    value -= (ax[later] - ax[i]).dot(ax[i] - ax[later])/(ax[later] - ax[i]).norm()/(ax[i] - ax[former]).norm();
        // }

        std::vector<Scalar> ax(GetRows());
        for (int i=0; i<GetRows(); ++i) ax[i] = x(i);
        auto fw = getADfun().Forward(0, ax);
        return fw[0];
    }

    void FillJacobianBlock(std::string var_name, Jacobian& jac) const override
    {
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        // std::vector<Vector2> ax(GetRows());
        // for (int i=0; i<horizon_; ++i) ax[i] = decorder_.decode(x(i), i);
        
        std::vector<Scalar> ax(GetRows());
        for (int i=0; i<GetRows(); ++i) ax[i] = x(i);
        auto ad_jac = getADfun().Jacobian(ax);
        for (int i=0; i<GetRows(); ++i)
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
            auto fx = decorder_.decodeX(ax[(i-1+horizon_)%horizon_], (i-1+horizon_)%horizon_);
            auto fy = decorder_.decodeY(ax[(i-1+horizon_)%horizon_], (i-1+horizon_)%horizon_);
            auto cx = decorder_.decodeX(ax[i], i);
            auto cy = decorder_.decodeY(ax[i], i);
            auto lx = decorder_.decodeX(ax[(i+1)%horizon_], (i+1)%horizon_);
            auto ly = decorder_.decodeY(ax[(i+1)%horizon_], (i+1)%horizon_);

            auto dot = (cx - fx)*(lx - cx) + (cy - fy)*(ly - cy);
            auto fnorm = CppAD::sqrt(CppAD::pow(cx - fx, 2.0) + CppAD::pow(cy - fy, 2.0));
            auto lnorm = CppAD::sqrt(CppAD::pow(lx - cx, 2.0) + CppAD::pow(ly - cy, 2.0));

            ay[0] = dot/fnorm/lnorm;
        }
        return CppAD::ADFun<Scalar>(ax, ay);
    }
    
    const int horizon_;
};
}
}