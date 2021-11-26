#ifndef SBS_MATH_ELASTICITY_H
#define SBS_MATH_ELASTICITY_H

#include "interpolation.h"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace sbs {
namespace math {
namespace differentiable {

template <class InterpolationFunctionType>
struct deformation_gradient_op_t
{
    using interpolate_op_type = InterpolationFunctionType;

    deformation_gradient_op_t(interpolate_op_type const& interpolate) : interpolate_op(interpolate)
    {
    }

    autodiff::Matrix3dual operator()(autodiff::Vector3dual X) const
    {
        using autodiff::at;
        using autodiff::gradient;
        using autodiff::jacobian;
        using autodiff::wrt;

        autodiff::Matrix3dual F;
        F.setZero();

        for (auto i = 0u; i < interpolate_op.uis.size(); ++i)
        {
            autodiff::Vector3dual gradphi = gradient(interpolate_op.phis[i], wrt(X), at(X));
            F += interpolate_op.uis[i] * gradphi.transpose();
        }

        return F;
    }
    autodiff::Matrix3dual operator()(autodiff::Vector3dual X, Eigen::Vector3dual& u) const
    {
        using autodiff::at;
        using autodiff::gradient;
        using autodiff::jacobian;
        using autodiff::wrt;

        u = interpolate_op(X);

        return (*this)(X);
    }

    interpolate_op_type const& interpolate_op;
};

template <class DeformationGradientFunctionType>
struct strain_op_t
{
    using deformation_gradient_op_type = DeformationGradientFunctionType;

    strain_op_t(deformation_gradient_op_type const& deformation_gradient)
        : deformation_gradient_op(deformation_gradient)
    {
    }

    autodiff::Matrix3dual
    operator()(autodiff::Vector3dual X, autodiff::Vector3dual& u, autodiff::Matrix3dual& F) const
    {
        F                       = deformation_gradient_op(X, u);
        autodiff::Matrix3dual E = (*this)(F);
        return E;
    }
    autodiff::Matrix3dual operator()(autodiff::Vector3dual X, autodiff::Matrix3dual& F) const
    {
        F                       = deformation_gradient_op(X);
        autodiff::Matrix3dual E = (*this)(F);
        return E;
    }

    autodiff::Matrix3dual operator()(autodiff::Matrix3dual F) const
    {
        autodiff::Matrix3dual I = autodiff::Matrix3dual::Identity();
        autodiff::Matrix3dual E = 0.5 * (F.transpose() * F - I);
        return E;
    }

    deformation_gradient_op_type const& deformation_gradient_op;
};

template <class StrainFunctionType>
struct stvk_strain_energy_density_op_t
{
    using strain_op_type = StrainFunctionType;

    stvk_strain_energy_density_op_t(
        strain_op_type const& strain,
        double young_modulus,
        double poisson_ratio)
        : strain_op(strain), mu(), lambda()
    {
        mu = (young_modulus) / (2. * (1. + poisson_ratio));
        lambda =
            (young_modulus * poisson_ratio) / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    }

    autodiff::dual operator()(
        autodiff::Vector3dual X,
        autodiff::Vector3dual& u,
        autodiff::Matrix3dual& F,
        autodiff::Matrix3dual& E) const
    {
        E = strain_op(X, u, F);
        return (*this)(E);
    }

    autodiff::dual
    operator()(autodiff::Vector3dual X, autodiff::Matrix3dual& F, autodiff::Matrix3dual& E) const
    {
        E = strain_op(X, F);
        return (*this)(E);
    }

    autodiff::dual operator()(autodiff::Matrix3dual F, autodiff::Matrix3dual& E) const
    {
        E = strain_op(F);
        return (*this)(E);
    }

    autodiff::dual operator()(autodiff::Matrix3dual E) const
    {
        autodiff::dual trace = E.trace();
        auto tr2             = trace * trace;
        autodiff::dual EdotE = (E.array() * E.array()).sum(); // contraction
        return mu * EdotE + 0.5 * lambda * tr2;
    }

    strain_op_type const& strain_op;
    double mu, lambda;
};

} // namespace differentiable

struct green_strain_op_t
{
    Eigen::Matrix3d operator()(Eigen::Matrix3d const& F) const
    {
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d E = 0.5 * (F.transpose() * F - I);
        return E;
    }
};

struct stvk_strain_energy_density_op_t
{
    stvk_strain_energy_density_op_t(scalar_type young_modulus, scalar_type poisson_ratio)
        : mu(), lambda()
    {
        mu = (young_modulus) / (2. * (1. + poisson_ratio));
        lambda =
            (young_modulus * poisson_ratio) / ((1. + poisson_ratio) * (1. - 2. * poisson_ratio));
    }

    scalar_type operator()(Eigen::Matrix3d const& E) const
    {
        scalar_type const trace = E.trace();
        auto tr2                = trace * trace;
        scalar_type const EdotE = (E.array() * E.array()).sum(); // contraction
        return mu * EdotE + 0.5 * lambda * tr2;
    }

    // Piola kirchhoff stress tensor: dPsi/dF
    Eigen::Matrix3d stress(Eigen::Matrix3d const& F, Eigen::Matrix3d const& E) const
    {
        scalar_type const tr        = E.trace();
        Eigen::Matrix3d const I     = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d const sigma = F * ((2. * mu * E) + (lambda * tr * I));
        return sigma;
    }

    scalar_type mu, lambda;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_ELASTICITY_H