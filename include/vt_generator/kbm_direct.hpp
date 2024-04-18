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
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Scalar> ax(x.size());
        for (int i=0; i<x.size(); ++i) ax[i] = x(i);

        VectorXd g(GetRows());
        for (int i=0; i<horizon_; ++i)
        {
            auto acc   = denormalizer_.denormalizeAcc(ax, i);
            auto steer = denormalizer_.denormalizeSteer(ax, i); 
            auto n     = denormalizer_.denormalizeN(ax, i);
            auto xi    = denormalizer_.denormalizeXi(ax, i);
            auto vx    = denormalizer_.denormalizeVx(ax, i);
            auto vy    = denormalizer_.denormalizeVy(ax, i);
            auto w     = denormalizer_.denormalizeW(ax, i);
            auto k     = curvature_.at(i);

            auto v    = std::sqrt(std::pow(vx,2.0) + std::pow(vy, 2.0));
            auto beta = std::atan(vy/vx);
            auto dtds = (1.0 - n*k) / (v*std::cos(xi + beta));
            auto ff   = -Cf_*((vy - Lf_*w)/vx - steer);
            auto fr   = -Cr_*(vy - Lr_*w)/vx;

            // state equation
            g(i*size_ + 0) = denormalizer_.denormalizeN(ax, (i+1)%horizon_)  - (n  + (std::tan(xi + beta)*(1-n*k))*ds_.at(i)); 
            g(i*size_ + 1) = denormalizer_.denormalizeXi(ax, (i+1)%horizon_) - (xi + (w - k*v*std::cos(xi + beta)/ (1.0 - n*k))*dtds*ds_.at(i));
            g(i*size_ + 2) = denormalizer_.denormalizeVx(ax, (i+1)%horizon_) - (vx + (acc   - ff*std::sin(steer)/m_ + vy*w)*dtds*ds_.at(i));
            g(i*size_ + 3) = denormalizer_.denormalizeVy(ax, (i+1)%horizon_) - (vy + (fr/m_ + ff*std::cos(steer)/m_ - vx*w)*dtds*ds_.at(i));
            g(i*size_ + 4) = denormalizer_.denormalizeW(ax, (i+1)%horizon_)  - (w  + (ff*Lf_*std::cos(steer) - fr*Lr_)*dtds*ds_.at(i));
            g(i*size_ + 5) = (ff + fr) / m_;
            g(i*size_ + 6) = beta;
        }

        return g;
    }

    void FillJacobianBlock(std::string var_name, Jacobian& jac) const override
    {
        auto x = GetVariables()->GetComponent("variables")->GetValues();
        std::vector<Scalar> ax(x.size());
        for (int i=0; i<x.size(); ++i) ax[i] = x(i);

        for (int i=0; i<horizon_; ++i)
        {
            auto acc   = denormalizer_.denormalizeAcc(ax, i);
            auto steer = denormalizer_.denormalizeSteer(ax, i); 
            auto n     = denormalizer_.denormalizeN(ax, i);
            auto xi    = denormalizer_.denormalizeXi(ax, i);
            auto vx    = denormalizer_.denormalizeVx(ax, i);
            auto vy    = denormalizer_.denormalizeVy(ax, i);
            auto w     = denormalizer_.denormalizeW(ax, i);
            auto k     = curvature_.at(i);
            auto ds    = ds_.at(i);
            auto beta  = std::atan(vy/vx);
            auto v     = std::sqrt(std::pow(vx,2.0) + std::pow(vy, 2.0));

            auto daccdnorm   = denormalizer_.daccdnorm(ax, i);
            auto dsteerdnorm = denormalizer_.dsteerdnorm(ax, i);
            auto dndnorm     = denormalizer_.dndnorm(ax, i);
            auto dxidnorm    = denormalizer_.dxidnorm(ax, i);
            auto dvxdnorm    = denormalizer_.dvxdnorm(ax, i);
            auto dvydnorm    = denormalizer_.dvydnorm(ax, i);
            auto dwdnorm     = denormalizer_.dwdnorm(ax, i);

            auto dbdvx =-vy/std::pow(vx, 2.0)/ (1 + std::pow(vy/vx, 2.0));
            auto dbdvy = 1.0/vx / (1 + std::pow(vy/vx, 2.0));
            auto dcndxi= -(1-n*k)*ds/std::pow(std::cos(xi + beta), 2.0);
            // dn/ds
            jac.coeffRef(i*size_ + 0, denormalizer_.n((i+1)%horizon_))  = denormalizer_.dndnorm(ax, (i+1)%horizon_);
            jac.coeffRef(i*size_ + 0, denormalizer_.n(i))  = (-1.0 + k*std::tan(xi + beta)*ds)*dndnorm; 
            jac.coeffRef(i*size_ + 0, denormalizer_.xi(i)) = dcndxi * dxidnorm; 
            jac.coeffRef(i*size_ + 0, denormalizer_.vx(i)) = dcndxi*dbdvx * dvxdnorm; 
            jac.coeffRef(i*size_ + 0, denormalizer_.vy(i)) = dcndxi*dbdvy * dvydnorm; 

            // ds/dt
            auto dsdt     = v *std::cos(beta + xi)   / (1-n*k);
            auto ds2dtdn  = v *std::cos(beta + xi)*k / std::pow(1-n*k, 2.0);
            auto ds2dtdxi =-v *std::sin(beta + xi)   / (1-n*k);
            auto ds2dtdvx = vx*std::cos(beta+xi)/(v*(1-n*k)) - v*std::sin(beta+xi)/(1-n*k) * dbdvx;
            auto ds2dtdvy = vy*std::cos(beta+xi)/(v*(1-n*k)) - v*std::sin(beta+xi)/(1-n*k) * dbdvy;

            // dxi
            jac.coeffRef(i*size_ + 1, denormalizer_.xi((i+1)%horizon_)) = denormalizer_.dxidnorm(ax, (i+1)%horizon_); 
            jac.coeffRef(i*size_ + 1, denormalizer_.n(i))  = -w*ds2dtdn*ds * dndnorm; 
            jac.coeffRef(i*size_ + 1, denormalizer_.xi(i)) =(-w*ds2dtdxi*ds - 1.0) * dxidnorm; 
            jac.coeffRef(i*size_ + 1, denormalizer_.vx(i)) = -w*ds2dtdvx*ds * dvxdnorm; 
            jac.coeffRef(i*size_ + 1, denormalizer_.vy(i)) = -w*ds2dtdvy*ds * dvydnorm; 
            jac.coeffRef(i*size_ + 1, denormalizer_.w(i))  = -dsdt*ds * dwdnorm; 
            
            // dFf/dx
            auto ff = -Cf_*((vy - Lf_*w)/vx - steer);
            auto dffdd  = Cf_;
            auto dffdvx = Cf_*(vy - Lf_*w)/std::pow(vx, 2.0);
            auto dffdvy =-Cf_/vx;
            auto dffdw  = Cf_*Lf_/vx;
            
            // d^2s/(dsdx)
            auto dtds     = (1-n*k)/(v*std::cos(xi + beta));
            auto dt2dsdn  = -k     /(v*std::cos(xi + beta));
            auto dt2dsdxi = (1-n*k)*std::sin(xi + beta)/(v*std::pow(std::cos(xi + beta), 2.0));
            auto dt2dsdvx =-(1-n*k)*std::sin(xi + beta)*vx/(std::pow(v, 3.0)*std::cos(xi + beta)) + dt2dsdxi*dbdvx;
            auto dt2dsdvy =-(1-n*k)*std::sin(xi + beta)*vy/(std::pow(v, 3.0)*std::cos(xi + beta)) + dt2dsdxi*dbdvy;

            // vx
            auto dvxdt = acc - ff*std::sin(steer)/m_ + vy*w;
            jac.coeffRef(i*size_ + 2, denormalizer_.vx((i+1)%horizon_)) = denormalizer_.dvxdnorm(ax, (i+1)%horizon_); 
            jac.coeffRef(i*size_ + 2, denormalizer_.acc(i))   =-dtds*ds * daccdnorm; 
            jac.coeffRef(i*size_ + 2, denormalizer_.steer(i)) = (std::sin(steer)*dffdd + ff*std::cos(steer))/m_*dtds*ds * dsteerdnorm; 
            jac.coeffRef(i*size_ + 2, denormalizer_.n(i))     =-dvxdt*dt2dsdn*ds * dndnorm; 
            jac.coeffRef(i*size_ + 2, denormalizer_.xi(i))    =-dvxdt*dt2dsdxi*ds * dxidnorm; 
            jac.coeffRef(i*size_ + 2, denormalizer_.vx(i))    =(-1.0 - (dvxdt*dt2dsdvx - dffdvx*std::sin(steer)/m_*dtds)*ds) * dvxdnorm; 
            jac.coeffRef(i*size_ + 2, denormalizer_.vy(i))    =-(dvxdt*dt2dsdvy - (dffdvy*std::sin(steer)/m_ - w)*dtds)*ds * dvydnorm; 
            jac.coeffRef(i*size_ + 2, denormalizer_.w(i))     = (dffdw*std::sin(steer)/m_ - vy)*dtds*ds * dwdnorm; 


            // fr
            auto fr = -Cr_*(vy - Lr_*w)/vx;
            auto dfrdvx = Cr_*(vy - Lr_*w)/std::pow(vx, 2.0);
            auto dfrdvy =-Cr_/vx;
            auto dfrdw  = Cr_*Lr_/vx;

            // dvy/ds
            auto dvydt = ((fr + ff*std::cos(steer))/m_ - vx*w);
            jac.coeffRef(i*size_ + 3, denormalizer_.vy((i+1)%horizon_)) = denormalizer_.dvydnorm(ax, (i+1)%horizon_); 
            jac.coeffRef(i*size_ + 3, denormalizer_.steer(i)) =-(dffdd*std::cos(steer) - ff*std::sin(steer))/m_*dtds*ds * dsteerdnorm; 
            jac.coeffRef(i*size_ + 3, denormalizer_.n(i))     =-dvydt*dt2dsdn*ds * dndnorm; 
            jac.coeffRef(i*size_ + 3, denormalizer_.xi(i))    =-dvydt*dt2dsdxi*ds * dxidnorm; 
            jac.coeffRef(i*size_ + 3, denormalizer_.vx(i))    = -(((dfrdvx + dffdvx*std::cos(steer))/m_ - w) *dtds + dvydt*dt2dsdvx)*ds * dvxdnorm; 
            jac.coeffRef(i*size_ + 3, denormalizer_.vy(i))    =(-(((dfrdvy + dffdvy*std::cos(steer))/m_)     *dtds + dvydt*dt2dsdvy)*ds - 1.0) * dvydnorm; 
            jac.coeffRef(i*size_ + 3, denormalizer_.w(i))     = -(((dfrdw  + dffdw *std::cos(steer))/m_ - vx)*dtds)*ds * dwdnorm; 

            // dw/ds
            auto dwdt = (ff*Lf_*std::cos(steer) - fr*Lr_)/Iz_;
            jac.coeffRef(i*size_ + 4, denormalizer_.w((i+1)%horizon_))  = denormalizer_.dwdnorm(ax, (i+1)%horizon_); 
            jac.coeffRef(i*size_ + 4, denormalizer_.steer(i)) =-(dffdd*Lf_*std::cos(steer) - ff*Lf_*std::sin(steer))/Iz_ * dtds*ds * dsteerdnorm; 
            jac.coeffRef(i*size_ + 4, denormalizer_.n(i))     =-dwdt*dt2dsdn*ds * dndnorm; 
            jac.coeffRef(i*size_ + 4, denormalizer_.xi(i))    =-dwdt*dt2dsdxi*ds * dxidnorm; 
            jac.coeffRef(i*size_ + 4, denormalizer_.vx(i))    = -((dffdvx*Lf_*std::cos(steer) - dfrdvx*Lr_)/Iz_*dtds + dwdt*dt2dsdvx)*ds * dvxdnorm; 
            jac.coeffRef(i*size_ + 4, denormalizer_.vy(i))    = -((dffdvy*Lf_*std::cos(steer) - dfrdvy*Lr_)/Iz_*dtds + dwdt*dt2dsdvy)*ds * dvydnorm; 
            jac.coeffRef(i*size_ + 4, denormalizer_.w(i))     =(-((dffdw* Lf_*std::cos(steer) - dfrdw *Lr_)/Iz_*dtds)*ds - 1.0) * dwdnorm;

            // lateral acc
            jac.coeffRef(i*size_ + 5, denormalizer_.steer(i)) = (dffdd)*dsteerdnorm;
            jac.coeffRef(i*size_ + 5, denormalizer_.vx(i))    = (dffdvx + dfrdvx)*dvxdnorm;
            jac.coeffRef(i*size_ + 5, denormalizer_.vy(i))    = (dffdvy + dfrdvy)*dvydnorm;
            jac.coeffRef(i*size_ + 5, denormalizer_.w(i))     = (dffdw + dfrdw)*dwdnorm;

            // slip angle
            jac.coeffRef(i*size_ + 6, denormalizer_.vx(i)) = (dbdvx)*dvxdnorm;
            jac.coeffRef(i*size_ + 6, denormalizer_.vy(i)) = (dbdvy)*dvydnorm;
            
        }
    }

private:
    const Scalar m_;
    const Scalar Iz_;
    const Scalar Cf_;
    const Scalar Cr_;
    const Scalar Lf_;
    const Scalar Lr_;

    const int size_;
    const int horizon_;
    const std::vector<Scalar> curvature_;
    const std::vector<Scalar> ds_;
    const VecBound bounds_;
};
}