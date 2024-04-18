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
        for (int i=0; i<x.size(); ++i) ax[i] = x(i);

        Scalar value = 0.0;
        for (int i=0; i<horizon_; ++i)
        {
            auto n  = denormalizer_.denormalizeN(ax, i);
            auto xi = denormalizer_.denormalizeXi(ax, i);
            auto vx = denormalizer_.denormalizeVx(ax, i);
            auto vy = denormalizer_.denormalizeVy(ax, i);
            auto& k  = curvature_[i];

            auto v = std::sqrt(std::pow(vx, 2) + std::pow(vy, 2));
            auto beta = std::atan(vy/vx);
            auto dtds = (1.0 - n*k) / (v*std::cos(xi + beta));
            // if (1.0 - n*k < 0) std::cout << "1.0 - n*k < 0: " << n << ", " << k << std::endl;
            // if (v < 0) std::cout << "v < 0: " << v << std::endl;
            // if (std::cos(xi + beta) < 0) std::cout << "cos(xi + beta) < 0: " << xi << ", " << beta << std::endl;
            value += dtds*ds_[i];
            value += std::pow(ax.at(denormalizer_.acc((i+1)%horizon_))   - ax.at(denormalizer_.acc(i)),   2.0);
            value += std::pow(ax.at(denormalizer_.steer((i+1)%horizon_)) - ax.at(denormalizer_.steer(i)), 2.0);
        }
        return value;
    }

    void FillJacobianBlock(std::string var_name, Jacobian& jac) const override
    {
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Scalar> ax(x.size());
        for (int i=0; i<x.size(); ++i)
        {
            ax[i] = x(i);
        }

        for (int i=0; i<horizon_; ++i)
        {
            auto n  = denormalizer_.denormalizeN(ax, i);
            auto xi = denormalizer_.denormalizeXi(ax, i);
            auto vx = denormalizer_.denormalizeVx(ax, i);
            auto vy = denormalizer_.denormalizeVy(ax, i);
            auto& k  = curvature_[i];
            auto& ds = ds_[i];

            auto v = std::sqrt(std::pow(vx, 2) + std::pow(vy, 2));
            auto beta = std::atan(vy/vx);
            auto dbdvx =-vy/std::pow(vx, 2.0)/ (1 + std::pow(vy/vx, 2.0));
            auto dbdvy = 1.0/vx / (1 + std::pow(vy/vx, 2.0));
            
            auto dtds     = (1-n*k)/(v*std::cos(xi + beta));
            auto dt2dsdn  = -k     /(v*std::cos(xi + beta));
            auto dt2dsdxi = (1-n*k)*std::sin(xi + beta)/(v*std::pow(std::cos(xi + beta), 2.0));
            auto dt2dsdvx =-(1-n*k)*std::sin(xi + beta)*vx/(std::pow(v, 3.0)*std::cos(xi + beta)) + dt2dsdxi*dbdvx;
            auto dt2dsdvy =-(1-n*k)*std::sin(xi + beta)*vy/(std::pow(v, 3.0)*std::cos(xi + beta)) + dt2dsdxi*dbdvy;

            jac.coeffRef(0, denormalizer_.acc(i))   = 4*ax.at(denormalizer_.acc(i))
                                                     -2*ax.at(denormalizer_.acc((i+1)%horizon_))
                                                     -2*ax.at(denormalizer_.acc((i-1+horizon_)%horizon_));
            jac.coeffRef(0, denormalizer_.steer(i)) = 4*ax.at(denormalizer_.steer(i))
                                                     -2*ax.at(denormalizer_.steer((i+1)%horizon_))
                                                     -2*ax.at(denormalizer_.steer((i-1+horizon_)%horizon_));
            jac.coeffRef(0, denormalizer_.n(i))     = dt2dsdn *ds*denormalizer_.dndnorm(ax, i);
            jac.coeffRef(0, denormalizer_.xi(i))    = dt2dsdxi*ds*denormalizer_.dxidnorm(ax, i);
            jac.coeffRef(0, denormalizer_.vx(i))    = dt2dsdvx*ds*denormalizer_.dvxdnorm(ax, i);
            jac.coeffRef(0, denormalizer_.vy(i))    = dt2dsdvy*ds*denormalizer_.dvydnorm(ax, i);
        }
    }

private:
    const int horizon_;
    std::vector<Scalar> ds_;
    std::vector<Scalar> curvature_;
};
}