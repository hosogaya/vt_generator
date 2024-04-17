#pragma once

#include <vt_generator/type.hpp>

namespace vt_generator
{

class KBM : public ifopt::ConstraintSet
{
public:
    KBM(const Scalar m, const Scalar Iz, 
        const Scalar Cf, const Scalar Cr,
        const Scalar Lf, const Scalar Lr, 
        const std::vector<Scalar>& curvature,
        const std::vector<Scalar>& ds,
        const VecBound& b, 
        const int horizon):
    ifopt::ConstraintSet(constraints_size_*horizon, "consraints_based_on_kbm"),
    m_(m), Iz_(Iz), Cf_(Cf), Cr_(Cr), Lf_(Lf), Lr_(Lr),
    curvature_(curvature), ds_(ds), bounds_(b), horizon_(horizon), size_(constraints_size_)
    {
    }
    ~KBM() {}

    VecBound GetBounds() const override
    {
        return bounds_;
    }

    VectorXd GetValues() const override
    {
        std::cout << "get value of constraints" << std::endl;
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Scalar> ax(x.size());
        for (int i=0; i<x.size(); ++i) ax[i] = x(i);
        std::vector<Scalar> fw = getADFun().Forward(0, ax);

        VectorXd g(GetRows());
        for (int i=0; i<GetRows(); ++i)
        {
            g(i) = fw[i];
        }
        return g;
    }

    void FillJacobianBlock(std::string var_name, Jacobian& jac) const override
    {
        std::cout << "fill jacobian of constraints" << std::endl;
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Scalar> ax(x.size());
        for (int i=0; i<x.size(); ++i) {
            ax[i] = x(i);
        }
        std::vector<Scalar> ad_jac = getADFun().Jacobian(ax);

        for (int i=0; i<horizon_; ++i)
        {
            // dn/ds
            jac.coeffRef(i*size_ + 0, denormalizer_.n(i))  = ad_jac[x.size()*(i*size_ + 0) + denormalizer_.n(i)]; 
            jac.coeffRef(i*size_ + 0, denormalizer_.xi(i)) = ad_jac[x.size()*(i*size_ + 0) + denormalizer_.xi(i)]; 
            jac.coeffRef(i*size_ + 0, denormalizer_.vx(i)) = ad_jac[x.size()*(i*size_ + 0) + denormalizer_.vx(i)]; 
            jac.coeffRef(i*size_ + 0, denormalizer_.vy(i)) = ad_jac[x.size()*(i*size_ + 0) + denormalizer_.vy(i)]; 

            // dxi/ds
            jac.coeffRef(i*size_ + 1, denormalizer_.n(i))  = ad_jac[x.size()*(i*size_ + 1) + denormalizer_.n(i)]; 
            jac.coeffRef(i*size_ + 1, denormalizer_.xi(i)) = ad_jac[x.size()*(i*size_ + 1) + denormalizer_.xi(i)]; 
            jac.coeffRef(i*size_ + 1, denormalizer_.vx(i)) = ad_jac[x.size()*(i*size_ + 1) + denormalizer_.vx(i)]; 
            jac.coeffRef(i*size_ + 1, denormalizer_.vy(i)) = ad_jac[x.size()*(i*size_ + 1) + denormalizer_.vy(i)]; 
            jac.coeffRef(i*size_ + 1, denormalizer_.w(i))  = ad_jac[x.size()*(i*size_ + 1) + denormalizer_.w(i)]; 
            
            // dvx/ds
            jac.coeffRef(i*size_ + 2, denormalizer_.acc(i))   = ad_jac[x.size()*(i*size_ + 2) + denormalizer_.acc(i)]; 
            jac.coeffRef(i*size_ + 2, denormalizer_.steer(i)) = ad_jac[x.size()*(i*size_ + 2) + denormalizer_.steer(i)]; 
            jac.coeffRef(i*size_ + 2, denormalizer_.n(i))     = ad_jac[x.size()*(i*size_ + 2) + denormalizer_.n(i)]; 
            jac.coeffRef(i*size_ + 2, denormalizer_.xi(i))    = ad_jac[x.size()*(i*size_ + 2) + denormalizer_.xi(i)]; 
            jac.coeffRef(i*size_ + 2, denormalizer_.vx(i))    = ad_jac[x.size()*(i*size_ + 2) + denormalizer_.vx(i)]; 
            jac.coeffRef(i*size_ + 2, denormalizer_.vy(i))    = ad_jac[x.size()*(i*size_ + 2) + denormalizer_.vy(i)]; 
            jac.coeffRef(i*size_ + 2, denormalizer_.w(i))     = ad_jac[x.size()*(i*size_ + 2) + denormalizer_.w(i)]; 


            // dvy/ds
            jac.coeffRef(i*size_ + 3, denormalizer_.steer(i)) = ad_jac[x.size()*(i*size_ + 3) + denormalizer_.steer(i)]; 
            jac.coeffRef(i*size_ + 3, denormalizer_.n(i))     = ad_jac[x.size()*(i*size_ + 3) + denormalizer_.n(i)]; 
            jac.coeffRef(i*size_ + 3, denormalizer_.xi(i))    = ad_jac[x.size()*(i*size_ + 3) + denormalizer_.xi(i)]; 
            jac.coeffRef(i*size_ + 3, denormalizer_.vx(i))    = ad_jac[x.size()*(i*size_ + 3) + denormalizer_.vx(i)]; 
            jac.coeffRef(i*size_ + 3, denormalizer_.vy(i))    = ad_jac[x.size()*(i*size_ + 3) + denormalizer_.vy(i)]; 
            jac.coeffRef(i*size_ + 3, denormalizer_.w(i))     = ad_jac[x.size()*(i*size_ + 3) + denormalizer_.w(i)]; 

            // dw/ds
            jac.coeffRef(i*size_ + 4, denormalizer_.steer(i)) = ad_jac[x.size()*(i*size_ + 4) + denormalizer_.steer(i)]; 
            jac.coeffRef(i*size_ + 4, denormalizer_.n(i))     = ad_jac[x.size()*(i*size_ + 4) + denormalizer_.n(i)]; 
            jac.coeffRef(i*size_ + 4, denormalizer_.xi(i))    = ad_jac[x.size()*(i*size_ + 4) + denormalizer_.xi(i)]; 
            jac.coeffRef(i*size_ + 4, denormalizer_.vx(i))    = ad_jac[x.size()*(i*size_ + 4) + denormalizer_.vx(i)]; 
            jac.coeffRef(i*size_ + 4, denormalizer_.vy(i))    = ad_jac[x.size()*(i*size_ + 4) + denormalizer_.vy(i)]; 
            jac.coeffRef(i*size_ + 4, denormalizer_.w(i))     = ad_jac[x.size()*(i*size_ + 4) + denormalizer_.w(i)]; 
            
            // i+1 step
            jac.coeffRef(i*size_ + 0, denormalizer_.n(i))  = ad_jac[x.size()*(i*size_ + 0) + denormalizer_.n((i+1)%horizon_)]; 
            jac.coeffRef(i*size_ + 1, denormalizer_.xi(i)) = ad_jac[x.size()*(i*size_ + 1) + denormalizer_.xi((i+1)%horizon_)]; 
            jac.coeffRef(i*size_ + 2, denormalizer_.vx(i)) = ad_jac[x.size()*(i*size_ + 2) + denormalizer_.vx((i+1)%horizon_)]; 
            jac.coeffRef(i*size_ + 3, denormalizer_.vy(i)) = ad_jac[x.size()*(i*size_ + 3) + denormalizer_.vy((i+1)%horizon_)]; 
            jac.coeffRef(i*size_ + 4, denormalizer_.w(i))  = ad_jac[x.size()*(i*size_ + 4) + denormalizer_.w((i+1)%horizon_)]; 


            // slip angle 
            jac.coeffRef(i*size_ + 5, denormalizer_.vx(i))    = ad_jac[x.size()*(i*size_ + 5) + denormalizer_.vx(i)]; 
            jac.coeffRef(i*size_ + 5, denormalizer_.vy(i))    = ad_jac[x.size()*(i*size_ + 5) + denormalizer_.vy(i)];

            // lateral acc
            jac.coeffRef(i*size_ + 6, denormalizer_.steer(i)) = ad_jac[x.size()*(i*size_ + 6) + denormalizer_.w(i)]; 
            jac.coeffRef(i*size_ + 6, denormalizer_.vx(i))    = ad_jac[x.size()*(i*size_ + 6) + denormalizer_.vx(i)]; 
            jac.coeffRef(i*size_ + 6, denormalizer_.vy(i))    = ad_jac[x.size()*(i*size_ + 6) + denormalizer_.vy(i)]; 
            jac.coeffRef(i*size_ + 6, denormalizer_.w(i))     = ad_jac[x.size()*(i*size_ + 6) + denormalizer_.w(i)]; 

            // input
            // jac.coeffRef(i*size_ + 7, denormalizer_.acc(i))     = ad_jac[x.size()*(i*size_ + 7) + denormalizer_.acc(i)]; 
            // jac.coeffRef(i*size_ + 7, denormalizer_.acc(i))     = ad_jac[x.size()*(i*size_ + 7) + denormalizer_.acc((i+1)%horizon_)]; 
            // jac.coeffRef(i*size_ + 8, denormalizer_.steer(i))   = ad_jac[x.size()*(i*size_ + 8) + denormalizer_.steer(i)]; 
            // jac.coeffRef(i*size_ + 8, denormalizer_.steer(i))   = ad_jac[x.size()*(i*size_ + 8) + denormalizer_.steer((i+1)%horizon_)]; 
        }
    }

private:
    CppAD::ADFun<Scalar> getADFun() const
    {
        // input: acc, steer
        // state: n, xi, vx, vy, w
        AdVector ay(size_*horizon_);
        AdVector ax(denormalizer_.size()*horizon_);
        CppAD::Independent(ax);

        // each step
        for (int i=0; i<horizon_; ++i)
        {
            // auto& acc   = ax.at(denormalizer_.acc(i));
            // auto& steer = ax.at(denormalizer_.steer(i));
            // auto& n     = ax.at(denormalizer_.n(i));
            // auto& xi    = ax.at(denormalizer_.xi(i));
            // auto& vx    = ax.at(denormalizer_.vx(i));
            // auto& vy    = ax.at(denormalizer_.vy(i));
            // auto& w     = ax.at(denormalizer_.w(i));
            // auto& k     = curvature_.at(i);

            auto acc   = denormalizer_.denormalizeAcc(ax, i);
            auto steer = denormalizer_.denormalizeSteer(ax, i); 
            auto n     = denormalizer_.denormalizeN(ax, i);
            auto xi    = denormalizer_.denormalizeXi(ax, i);
            auto vx    = denormalizer_.denormalizeVx(ax, i);
            auto vy    = denormalizer_.denormalizeVy(ax, i);
            auto w     = denormalizer_.denormalizeW(ax, i);
            auto k     = curvature_.at(i);
            // using soft max
            auto vx_soft = 1.0/(vx + CppAD::log(1 + CppAD::exp(-2*vx*alpha_))/alpha_);
            auto v = CppAD::pow(vx, 2) + CppAD::pow(vy, 2);
            auto beta = CppAD::atan2(vy, vx_soft);
            auto ds = (v*CppAD::cos(xi + beta));
            auto dtds = (1.0 - n*k) / (ds + CppAD::log(1 + CppAD::exp(-2*ds*alpha_))/alpha_);
            auto ff = -Cf_*((vy - Lf_*w)/vx_soft - steer);
            auto fr = -Cr_*(vy - Lr_*w)/vx_soft;

            // auto v = CppAD::pow(vx, 2) + CppAD::pow(vy, 2);
            // auto beta = CppAD::atan2(vy, vx);
            // auto dtds = (1.0 + n*k) / (v*CppAD::cos(xi + beta));
            // auto ff = -Cf_*((vy - Lf_*w)/vx - steer);
            // auto fr = -Cr_*(vy - Lr_*w)/vx;

            // std::cout << "ff: " << ff << std::endl;
            // std::cout << "fr: " << fr << std::endl;
            // std::cout << "dtds: " << dtds << std::endl;
            // std::cout << "beta: " << beta << std::endl;
            // std::cout << "v: " << v << std::endl;
            // std::cout << "vx: " << vx << std::endl;

            // state equation
            ay[i*size_ + 0] = denormalizer_.denormalizeN(ax, (i+1)%horizon_) - (n  + (CppAD::tan(xi + beta)*(1-n*k))*ds_.at(i)); 
            ay[i*size_ + 1] = denormalizer_.denormalizeXi(ax, (i+1)%horizon_) - (xi + (w - k*v*CppAD::cos(xi + beta)/ (1.0 - n*k))*dtds*ds_.at(i));
            ay[i*size_ + 2] = denormalizer_.denormalizeVx(ax, (i+1)%horizon_) - (vx + (acc - ff*CppAD::sin(steer)/m_ + vy*w)*dtds*ds_.at(i));
            ay[i*size_ + 3] = denormalizer_.denormalizeVy(ax, (i+1)%horizon_) - (vy + (fr/m_ + ff*CppAD::cos(steer)/m_ - vx*w)*dtds*ds_.at(i));
            ay[i*size_ + 4] = denormalizer_.denormalizeW(ax, (i+1)%horizon_) - (w  + (ff*Lf_*CppAD::cos(steer) - fr*Lr_)*dtds*ds_.at(i));
            
            // slip angle
            ay[i*size_ + 5] = beta;
            // lateral acc
            ay[i*size_ + 6] = (ff + fr)/m_;

            // input
            // ay[i*size_ + 7] = denormalizer_.denormalizeAcc(ax, (i+1)%horizon_) - acc;
            // ay[i*size_ + 8] = denormalizer_.denormalizeSteer(ax, (i+1)%horizon_) - steer;

            // for (int j=0; j<size_; ++j)
            // {
            //     if (CppAD::isnan(ay[i*size_+j]))
            //     {
            //         std::cout << "(" << i << ", " << j << ")" << ", ";
            //     }
            // }
        }
        CppAD::ADFun<Scalar> f(ax, ay);
        return f;
    }

private:
    const Scalar m_;
    const Scalar Iz_;
    const Scalar Cf_;
    const Scalar Cr_;
    const Scalar Lf_;
    const Scalar Lr_;

    const int size_ = 7;
    const int horizon_;
    const std::vector<Scalar> curvature_;
    const std::vector<Scalar> ds_;
    const VecBound bounds_;
    const Scalar alpha_ = 1e3;
};
}