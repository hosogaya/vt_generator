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
    ifopt::ConstraintSet(9*horizon, "consraints_based_on_kbm"),
    m_(m), Iz_(Iz), Cf_(Cf), Cr_(Cr), Lf_(Lf), Lr_(Lr),
    curvature_(curvature), ds_(ds), bounds_(b), horizon_(horizon)
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
        for (int i=0; i<x.size(); ++i) {
            ax[i] = x(i);
        }
        std::vector<Scalar> fw = getADFun(x).Forward(0, ax);

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
        std::vector<Scalar> ad_jac = getADFun(x).Jacobian(ax);

        for (int i=0; i<horizon_; ++i)
        {
            // dn/ds
            jac.coeffRef(i*size_ + 0, var_index.n(i))  = ad_jac[x.size()*(i*size_ + 0) + var_index.n(i)]; 
            jac.coeffRef(i*size_ + 0, var_index.xi(i)) = ad_jac[x.size()*(i*size_ + 0) + var_index.xi(i)]; 
            jac.coeffRef(i*size_ + 0, var_index.vx(i)) = ad_jac[x.size()*(i*size_ + 0) + var_index.vx(i)]; 
            jac.coeffRef(i*size_ + 0, var_index.vy(i)) = ad_jac[x.size()*(i*size_ + 0) + var_index.vy(i)]; 

            // dxi/ds
            jac.coeffRef(i*size_ + 1, var_index.n(i))  = ad_jac[x.size()*(i*size_ + 1) + var_index.n(i)]; 
            jac.coeffRef(i*size_ + 1, var_index.xi(i)) = ad_jac[x.size()*(i*size_ + 1) + var_index.xi(i)]; 
            jac.coeffRef(i*size_ + 1, var_index.vx(i)) = ad_jac[x.size()*(i*size_ + 1) + var_index.vx(i)]; 
            jac.coeffRef(i*size_ + 1, var_index.vy(i)) = ad_jac[x.size()*(i*size_ + 1) + var_index.vy(i)]; 
            jac.coeffRef(i*size_ + 1, var_index.w(i))  = ad_jac[x.size()*(i*size_ + 1) + var_index.w(i)]; 
            
            // dvx/ds
            jac.coeffRef(i*size_ + 2, var_index.acc(i))   = ad_jac[x.size()*(i*size_ + 2) + var_index.acc(i)]; 
            jac.coeffRef(i*size_ + 2, var_index.steer(i)) = ad_jac[x.size()*(i*size_ + 2) + var_index.steer(i)]; 
            jac.coeffRef(i*size_ + 2, var_index.n(i))     = ad_jac[x.size()*(i*size_ + 2) + var_index.n(i)]; 
            jac.coeffRef(i*size_ + 2, var_index.xi(i))    = ad_jac[x.size()*(i*size_ + 2) + var_index.xi(i)]; 
            jac.coeffRef(i*size_ + 2, var_index.vx(i))    = ad_jac[x.size()*(i*size_ + 2) + var_index.vx(i)]; 
            jac.coeffRef(i*size_ + 2, var_index.vy(i))    = ad_jac[x.size()*(i*size_ + 2) + var_index.vy(i)]; 
            jac.coeffRef(i*size_ + 2, var_index.w(i))     = ad_jac[x.size()*(i*size_ + 2) + var_index.w(i)]; 


            // dvy/ds
            jac.coeffRef(i*size_ + 3, var_index.steer(i)) = ad_jac[x.size()*(i*size_ + 3) + var_index.steer(i)]; 
            jac.coeffRef(i*size_ + 3, var_index.n(i))     = ad_jac[x.size()*(i*size_ + 3) + var_index.n(i)]; 
            jac.coeffRef(i*size_ + 3, var_index.xi(i))    = ad_jac[x.size()*(i*size_ + 3) + var_index.xi(i)]; 
            jac.coeffRef(i*size_ + 3, var_index.vx(i))    = ad_jac[x.size()*(i*size_ + 3) + var_index.vx(i)]; 
            jac.coeffRef(i*size_ + 3, var_index.vy(i))    = ad_jac[x.size()*(i*size_ + 3) + var_index.vy(i)]; 
            jac.coeffRef(i*size_ + 3, var_index.w(i))     = ad_jac[x.size()*(i*size_ + 3) + var_index.w(i)]; 

            // dw/ds
            jac.coeffRef(i*size_ + 4, var_index.steer(i)) = ad_jac[x.size()*(i*size_ + 4) + var_index.steer(i)]; 
            jac.coeffRef(i*size_ + 4, var_index.n(i))     = ad_jac[x.size()*(i*size_ + 4) + var_index.n(i)]; 
            jac.coeffRef(i*size_ + 4, var_index.xi(i))    = ad_jac[x.size()*(i*size_ + 4) + var_index.xi(i)]; 
            jac.coeffRef(i*size_ + 4, var_index.vx(i))    = ad_jac[x.size()*(i*size_ + 4) + var_index.vx(i)]; 
            jac.coeffRef(i*size_ + 4, var_index.vy(i))    = ad_jac[x.size()*(i*size_ + 4) + var_index.vy(i)]; 
            jac.coeffRef(i*size_ + 4, var_index.w(i))     = ad_jac[x.size()*(i*size_ + 4) + var_index.w(i)]; 
            
            // i+1 step
            jac.coeffRef(i*size_ + 0, var_index.n(i))  = ad_jac[x.size()*(i*size_ + 0) + var_index.n((i+1)%horizon_)]; 
            jac.coeffRef(i*size_ + 1, var_index.xi(i)) = ad_jac[x.size()*(i*size_ + 1) + var_index.xi((i+1)%horizon_)]; 
            jac.coeffRef(i*size_ + 2, var_index.vx(i)) = ad_jac[x.size()*(i*size_ + 2) + var_index.vx((i+1)%horizon_)]; 
            jac.coeffRef(i*size_ + 3, var_index.vy(i)) = ad_jac[x.size()*(i*size_ + 3) + var_index.vy((i+1)%horizon_)]; 
            jac.coeffRef(i*size_ + 4, var_index.w(i))  = ad_jac[x.size()*(i*size_ + 4) + var_index.w((i+1)%horizon_)]; 


            // slip angle 
            jac.coeffRef(i*size_ + 5, var_index.vx(i))    = ad_jac[x.size()*(i*size_ + 5) + var_index.vx(i)]; 
            jac.coeffRef(i*size_ + 5, var_index.vy(i))    = ad_jac[x.size()*(i*size_ + 5) + var_index.vy(i)];

            // lateral acc
            jac.coeffRef(i*size_ + 6, var_index.steer(i)) = ad_jac[x.size()*(i*size_ + 6) + var_index.w(i)]; 
            jac.coeffRef(i*size_ + 6, var_index.vx(i))    = ad_jac[x.size()*(i*size_ + 6) + var_index.vx(i)]; 
            jac.coeffRef(i*size_ + 6, var_index.vy(i))    = ad_jac[x.size()*(i*size_ + 6) + var_index.vy(i)]; 
            jac.coeffRef(i*size_ + 6, var_index.w(i))     = ad_jac[x.size()*(i*size_ + 6) + var_index.w(i)]; 

            // input
            jac.coeffRef(i*size_ + 7, var_index.acc(i))     = ad_jac[x.size()*(i*size_ + 7) + var_index.acc(i)]; 
            jac.coeffRef(i*size_ + 7, var_index.acc(i))     = ad_jac[x.size()*(i*size_ + 7) + var_index.acc((i+1)%horizon_)]; 
            jac.coeffRef(i*size_ + 8, var_index.steer(i))   = ad_jac[x.size()*(i*size_ + 8) + var_index.steer(i)]; 
            jac.coeffRef(i*size_ + 8, var_index.steer(i))   = ad_jac[x.size()*(i*size_ + 8) + var_index.steer((i+1)%horizon_)]; 
        }
    }

private:
    CppAD::ADFun<Scalar> getADFun(const VectorXd& x) const
    {
        // input: acc, steer
        // state: n, xi, vx, vy, w
        AdVector ay(size_*horizon_);
        AdVector ax(x.size());
        CppAD::Independent(ax);

        // each step
        for (int i=0; i<horizon_; ++i)
        {
            auto& acc   = ax.at(var_index.acc(i));
            auto& steer = ax.at(var_index.steer(i));
            auto& n     = ax.at(var_index.n(i));
            auto& xi    = ax.at(var_index.xi(i));
            auto& vx    = ax.at(var_index.vx(i));
            auto& vy    = ax.at(var_index.vy(i));
            auto& w     = ax.at(var_index.w(i));
            auto& k     = curvature_.at(i);
            
            // using soft max
            auto vx_inv = 1.0/(vx + CppAD::log(1 + CppAD::exp(-2*vx*alpha_))/alpha_);
            auto v = CppAD::pow(vx, 2) + CppAD::pow(vy, 2);
            auto beta = CppAD::atan(vy*vx_inv);
            auto ds = (v*CppAD::cos(xi + beta));
            auto dtds = (1.0 - n*k) / (ds + CppAD::log(1 + CppAD::exp(-2*ds*alpha_))/alpha_);
            auto ff = -Cf_*((vy - Lf_*w)*vx_inv - steer);
            auto fr = -Cr_*(vy - Lr_*w)*vx_inv;

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
            ay[i*size_ + 0] = ax.at(var_index.n((i+1)%horizon_))  - (n + (v*CppAD::sin(xi + beta))*dtds*ds_.at(i)); 
            ay[i*size_ + 1] = ax.at(var_index.xi((i+1)%horizon_)) - (xi + (w - k*v*CppAD::cos(xi + beta)/ (1.0 - n*k))*dtds*ds_.at(i));
            ay[i*size_ + 2] = ax.at(var_index.vx((i+1)%horizon_)) - (vx + (acc - ff*CppAD::sin(steer)/m_ + vy*w)*dtds*ds_.at(i));
            ay[i*size_ + 3] = ax.at(var_index.vy((i+1)%horizon_)) - (vy + (fr/m_ + ff*CppAD::cos(steer)/m_ - vx*w)*dtds*ds_.at(i));
            ay[i*size_ + 4] = ax.at(var_index.w((i+1)%horizon_))  - (w + (ff*Lf_*CppAD::cos(steer) - fr*Lr_)*dtds*ds_.at(i));
            
            // slip angle
            ay[i*size_ + 5] = beta;
            // lateral acc
            ay[i*size_ + 6] = (ff + fr)/m_;

            // input
            ay[i*size_ + 7] = ax[var_index.acc((i+1)%horizon_)] - acc;
            ay[i*size_ + 8] = ax[var_index.steer((i+1)%horizon_)] - steer;

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

    const int size_ = 9;
    const int horizon_;
    const std::vector<Scalar> curvature_;
    const std::vector<Scalar> ds_;
    const VecBound bounds_;
    const Scalar alpha_ = 1e3;
};
}