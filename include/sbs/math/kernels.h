#ifndef SBS_MATH_KERNELS_H
#define SBS_MATH_KERNELS_H

#include "sbs/aliases.h"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace sbs {
namespace math {
namespace differentiable {

struct poly6_kernel_t
{
    poly6_kernel_t() = default;
    poly6_kernel_t(autodiff::Vector3dual Xi, scalar_type h)
        : Xi(Xi), h(h), h2(h * h), h9(h * h * h * h * h * h * h * h * h), alpha()
    {
        alpha = 315. / (64. * pi * h9);
    }

    autodiff::dual operator()(autodiff::Vector3dual X) const
    {
        autodiff::dual r    = (X - Xi).norm();
        autodiff::dual r2   = r * r;
        autodiff::dual h_r  = (h2 - r2);
        autodiff::dual h_r3 = h_r * h_r * h_r;
        return alpha * h_r3;
    }

    Eigen::Vector3d grad(Eigen::Vector3d const& X) const
    {
        Eigen::Vector3d const diff           = X - Xi.cast<scalar_type>();
        scalar_type const r                  = diff.norm();
        scalar_type const r2                 = r * r;
        scalar_type constexpr mult_315_64_pi = 315. / (64. * 3.14159265359);
        scalar_type const mult = -6. * r * mult_315_64_pi / static_cast<scalar_type>(h9);
        scalar_type const hr3  = (static_cast<scalar_type>(h2) - r2) *
                                (static_cast<scalar_type>(h2) - r2) *
                                (static_cast<scalar_type>(h2) - r2);
        scalar_type const dWdr      = mult * hr3;
        Eigen::Vector3d const drdxj = diff.normalized();
        return dWdr * drdxj;
    }

    static constexpr scalar_type pi = 3.1415926535;

    autodiff::Vector3dual Xi;
    autodiff::dual h;
    autodiff::dual h2;
    autodiff::dual h9;
    autodiff::dual alpha;
};

} // namespace differentiable

struct poly6_kernel_t
{
    poly6_kernel_t() = default;
    poly6_kernel_t(Eigen::Vector3d const& Xi, scalar_type h)
        : Xi(Xi), h(h), h2(h * h), h9(h * h * h * h * h * h * h * h * h), alpha()
    {
        alpha = 315. / (64. * pi * h9);
    }

    scalar_type operator()(Eigen::Vector3d const& X) const { return eval(X); }
    scalar_type eval(Eigen::Vector3d const& X) const
    {
        scalar_type const r    = (X - Xi).norm();
        scalar_type const r2   = r * r;
        scalar_type const h_r  = (h2 - r2);
        scalar_type const h_r3 = h_r * h_r * h_r;
        return alpha * h_r3;
    }

    Eigen::Vector3d grad(Eigen::Vector3d const& X) const
    {
        Eigen::Vector3d const diff           = X - Xi;
        scalar_type const r                  = diff.norm();
        scalar_type const r2                 = r * r;
        scalar_type constexpr mult_315_64_pi = 315. / (64. * pi);
        scalar_type const mult               = -6. * r * mult_315_64_pi / h9;
        scalar_type const hr3                = (h2 - r2) * (h2 - r2) * (h2 - r2);
        scalar_type const dWdr               = mult * hr3;
        Eigen::Vector3d const drdxj          = diff.normalized();
        return dWdr * drdxj;
    }

    Eigen::Matrix3d hessian(Eigen::Vector3d const& X) const
    {
        Eigen::Vector3d const diff      = X - Xi;
        scalar_type const r             = diff.norm();
        scalar_type const r2            = r * r;
        scalar_type const r3            = r2 * r;
        Eigen::Matrix3d const I         = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d const c1        = (1. / r) * I;
        Eigen::Matrix3d const c2        = (1. / r3) * (diff * diff.transpose());
        Eigen::Matrix3d const hessian_r = c1 - c2;

        Eigen::Vector3d const dr_dX          = diff / r;
        scalar_type constexpr mult_315_64_pi = 315. / (64. * 3.14159265359);
        scalar_type const mult               = -6. * r * mult_315_64_pi / h9;
        scalar_type const hr3                = (h2 - r2) * (h2 - r2) * (h2 - r2);
        scalar_type const dW_dr              = mult * hr3;

        scalar_type const d2W_dr2_lhs = (24. * mult_315_64_pi * r2 * (h2 - r2)) / h9;
        scalar_type const d2W_dr2_rhs = (6. * mult_315_64_pi * (h2 - r2) * (h2 - r2)) / h9;
        scalar_type const d2W_dr2     = d2W_dr2_lhs - d2W_dr2_rhs;

        Eigen::Matrix3d H{};
        H.setZero();
        for (int m = 0; m < 3; ++m)
        {
            for (int n = 0; n < 3; ++n)
            {
                H(m, n) = (dW_dr * hessian_r(m, n)) + (d2W_dr2 * dr_dX(n) * dr_dX(n));
            }
        }

        return H;
    }

    static constexpr scalar_type pi = 3.1415926535;

    Eigen::Vector3d Xi;
    scalar_type h;
    scalar_type h2;
    scalar_type h9;
    scalar_type alpha;
};

struct quartic_spline_kernel_t
{
    quartic_spline_kernel_t() = default;
    quartic_spline_kernel_t(autodiff::Vector3dual Xi, scalar_type r) : Xi(Xi), r(r) {}

    autodiff::dual operator()(autodiff::Vector3dual X) const
    {
        autodiff::dual d = (X - Xi).norm() / r;
        if (d > 1)
            return 0.;

        autodiff::dual d2 = d * d;
        autodiff::dual d3 = d2 * d;
        autodiff::dual d4 = d2 * d2;
        return 1. - (6. * d2) + (8. * d3) - (3. * d4);
    }

    autodiff::Vector3dual Xi;
    autodiff::dual r;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_KERNELS_H