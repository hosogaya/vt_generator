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
        Scalar value = 0.0;
        std::vector<Vector2> ax(x.size());
        for (int i=0; i<horizon_; ++i) ax[i] = decorder_.decode(x(i), i);
        for (int i=0; i<horizon_; ++i)
        {
           int former = (i-1+horizon_)%horizon_;
           int later  = (i+1)%horizon_;
           value -= (ax[later] - ax[i]).dot(ax[i] - ax[later])/(ax[later] - ax[i]).norm()/(ax[i] - ax[former]).norm();
        }

        // std::vector<Scalar> ax(x.size());
        // for (int i=0; i<x.size(); ++i) ax[i] = x(i);
        // auto fw = getADfun().Forward(0, ax);
        return value;
    }

    void FillJacobianBlock(std::string var_name, Jacobian& jac) const override
    {
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Vector2> ax(x.size());
        for (int i=0; i<horizon_; ++i) ax[i] = decorder_.decode(x(i), i);
        
        for (int i=0; i<horizon_; ++i)
        {
            int kp1 = (i+1)%horizon_;
            int kp2 = (i+2)%horizon_;
            int km1 = (i-1+horizon_)%horizon_;
            int km2 = (i-2+horizon_)%horizon_;

            auto xp1 = ax[kp1];
            auto xp2 = ax[kp2];
            auto x   = ax[i];
            auto xm1 = ax[km1];
            auto xm2 = ax[km2];
            
            auto np1 = (xp1 - x).dot(xp2 - xp1);
            auto n   = (x - xm1).dot(xp1 - x);
            auto nm1 = (xm1 - xm2).dot(x - xm1);

            auto d1p1 = (xp1 - x).norm();
            auto d2p1 = (xp2 - xp1).norm();
            auto d1   = (x - xm1).norm();
            auto d2   = (xp1 - x).norm();
            auto d1m1 = (xm1 - xm2).norm();
            auto d2m1 = (x - xm1).norm();

            auto ap1  = (xp1 - xp2)/(d1p1*d2p1)   + (xp1 - x)/std::pow(d1p1, 3) * np1/d2p1;
            auto a    = (xp1 - 2*x + xm1)/(d1*d2) + (x - xp1)/std::pow(d1  , 3) * n  /d2    + (x - xm1)/std::pow(d2, 3)*n/d1;
            auto am1  = (xm1 - xm2)/(d1m1*d2m1)   + (xm1 - x)/std::pow(d2m1, 3) * nm1/d1m1;

            auto sum = ap1 + a + am1;
            auto dxdn = decorder_.dxdn(i);
            auto dydn = decorder_.dydn(i);
            jac.coeffRef(0, i) = -sum.x()*dxdn - sum.y()*dydn;
        }

        // std::vector<Scalar> ax(x.size());
        // for (int i=0; i<x.size(); ++i) ax[i] = x(i);
        // auto ad_jac = getADfun().Jacobian(ax);
        // for (int i=0; i<x.size(); ++i)
        // {
        //     jac.coeffRef(0, i) = ad_jac.at(i);
        // }
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
            
            ay[0] -= dot/fnorm/lnorm;
        }
        return CppAD::ADFun<Scalar>(ax, ay);
    }
    
    const int horizon_;
};
}
}