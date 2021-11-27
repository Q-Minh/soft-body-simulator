#ifndef SBS_MATH_MLS_H
#define SBS_MATH_MLS_H

#include "basis_functions.h"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <vector>

namespace sbs {
namespace math {

template <class KernelFunctionType, unsigned int Order>
struct mls_interpolation_op_t
{
    using kernel_function_type = KernelFunctionType;
    using self_type            = mls_interpolation_op_t<KernelFunctionType, Order>;
    using basis_function_type  = differentiable::mls_basis_function_t<kernel_function_type, Order>;

    mls_interpolation_op_t() = default;
    mls_interpolation_op_t(
        std::vector<autodiff::Vector3dual> const& uis,
        std::vector<basis_function_type> const& phis)
        : uis(uis), phis(phis), gradphis()
    {
    }

    autodiff::Vector3dual operator()(autodiff::Vector3dual X) const
    {
        autodiff::Vector3dual u{0., 0., 0.};
        for (auto i = 0u; i < uis.size(); ++i)
        {
            basis_function_type const& phi = phis[i];
            u += uis[i] * phi(X);
        }
        return u;
    }

    void cache_grad_phis(autodiff::Vector3dual X)
    {
        // Cache basis function gradients
        gradphis.clear();
        gradphis.reserve(uis.size());
        for (auto i = 0u; i < uis.size(); ++i)
        {
            using autodiff::at;
            using autodiff::gradient;
            using autodiff::wrt;

            basis_function_type const& phi = phis[i];
            autodiff::Vector3dual gradphi  = gradient(phi, wrt(X), at(X));
            gradphis.push_back(gradphi);
        }
    }

    std::vector<autodiff::Vector3dual> uis;      ///< Interpolation coefficients
    std::vector<basis_function_type> phis;       ///< MLS basis functions
    std::vector<autodiff::Vector3dual> gradphis; ///< Cached basis function gradients
};

template <class MlsInterpolationType>
struct mls_deformation_gradient_op_t
{
    using mls_interpolation_op_type = MlsInterpolationType;

    mls_deformation_gradient_op_t() = default;
    mls_deformation_gradient_op_t(mls_interpolation_op_type const& op) : interpolate_op(op) {}

    /**
     * @brief
     * Computes the deformation gradient at X using an MLS interpolation.
     * Assumes the basis functiong gradients (gradphis) have already been precomputed on the
     * mls_interpolation_op_type.
     * @param X
     * @return
     */
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
            autodiff::Vector3dual const& gradphi = interpolate_op.gradphis[i];
            F += interpolate_op.uis[i] * gradphi.transpose();
        }

        return F;
    }

    mls_interpolation_op_type const& interpolate_op;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_MLS_H