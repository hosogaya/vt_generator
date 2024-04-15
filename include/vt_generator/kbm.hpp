#pragma once

#include <vt_generator/type.hpp>

namespace vt_generator
{
class KBM : ifopt::ConstraintSet
{
public:
    KBM(const Scalar m, const Scalar Iz, 
        const Scalar Cf, const Scalar Cr,
        const Scalar Lf, const Scalar Lr, 
        const std::vector<Scalar>& curvature,
        const VecBound& b, 
        const int step_size):
    ifopt::ConstraintSet(1, "consraints_based_on_kbm"),
    m_(m), Iz_(Iz), Cf_(Cf), Cr_(Cr), Lf_(Lf), Lr_(Lr), curvature_(curvature), step_size_(step_size)
    {
        bounds_.resize(size_*step_size_);
        for (int i=0; i<step_size_; ++i)
        {
            for (int j=0; j<size_; ++j)
            {
                bounds_[i*size_ + j] = b[j];
            }
        }
    }
    ~KBM() {}

    VecBound GetBounds() const override
    {
        return bounds_;
    }

    VectorXd GetValues() const override
    {
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Scalar> ax(x.size());
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
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Scalar> ax(x.size());
        std::vector<Scalar> ad_jac = getADFun().Jacobian(ax);
        for (int i=0; i<GetRows(); ++i)
        {
            for (int j=0; j<x.size(); ++j)
            {
                jac.coeffRef(i,j) = ad_jac[i*x.size() + j];
            }
        }
    }


    CppAD::ADFun<Scalar> getADFun() const
    {
        // input: acc, delta
        // state: n, xi, vx, vy, w
        AdVector ay(size_*step_size_);
        AdVector ax(ind_.size()*step_size_);
        CppAD::Independent(ax);
        
        // each step
        for (int i=0; i<step_size_; ++i)
        {
            auto& acc = ax[ind_.acc(i)];
            auto& delta = ax[ind_.delta(i)];
            auto& n = ax[ind_.n(i)];
            auto& xi = ax[ind_.xi(i)];
            auto& vx = ax[ind_.vx(i)];
            auto& vy = ax[ind_.vy(i)];
            auto& w = ax[ind_.w(i)];
            auto& k = curvature_[i];

            auto v = CppAD::pow(vx, 2) + CppAD::pow(vy, 2);
            auto beta = CppAD::atan2(vy, vx);
            auto dtds = (1.0 + n*k) / v*CppAD::cos(xi + beta);
            auto ff = -Cf_*((vy - Lf_*w)/vx - delta);
            auto fr = -Cr_*(vy - Lr_*w)/vx;

            // state equation
            ay[i*size_ + 0] = ax[ind_.n((i+1)%step_size_)]  - (n + (v*CppAD::sin(xi + beta))*dtds); 
            ay[i*size_ + 1] = ax[ind_.xi((i+1)%step_size_)] - (xi + (w - k*v*CppAD::cos(xi + beta)/ (1.0 + n*k))*dtds);
            ay[i*size_ + 2] = ax[ind_.vx((i+1)%step_size_)] - (vx + (acc - ff*CppAD::sin(delta)/m_ + vy*w)*dtds);
            ay[i*size_ + 3] = ax[ind_.vy((i+1)%step_size_)] - (vy + (fr/m_ + ff*CppAD::cos(delta)/m_ - vx*w)*dtds);
            ay[i*size_ + 4] = ax[ind_.w((i+1)%step_size_)]  - (w + (ff*Lf_*CppAD::cos(delta) - fr*Lr_)*dtds);
            
            // slip angle
            ay[i*size_ + 5] = beta;
            // lateral acc
            ay[i*size_ + 6] = (ff + fr)/m_;

            // input
            ay[i*size_ + 7] = ax[ind_.acc((i+1)%step_size_)] - acc;
            ay[i*size_ + 8] = ax[ind_.delta((i+1)%step_size_)] - delta;
        }
        return CppAD::ADFun<Scalar>(ax, ay);
    }

private:
    const Scalar m_;
    const Scalar Iz_;
    const Scalar Cf_;
    const Scalar Cr_;
    const Scalar Lf_;
    const Scalar Lr_;

    Index ind_;
    const int size_ = 10;
    const int step_size_;
    std::vector<Scalar> curvature_;
    VecBound bounds_;
};
}