#ifndef SBS_MATH_INTERPOLATION_H
#define SBS_MATH_INTERPOLATION_H

#include <array>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <vector>

namespace sbs {
namespace math {

template <class BasisFunctionType, int NumBasisFunctions = -1>
struct interpolation_op_t
{
    using basis_function_op_type = BasisFunctionType;

    interpolation_op_t(
        std::array<autodiff::Vector3dual, NumBasisFunctions> const& uis,
        std::array<basis_function_op_type, NumBasisFunctions> const& phis)
        : uis(uis), phis(phis)
    {
    }

    autodiff::Vector3dual operator()(autodiff::Vector3dual const& X) const
    {
        autodiff::Vector3dual u;
        u << 0., 0., 0.;

        for (auto i = 0u; i < uis.size(); ++i)
            u += uis[i] * phis[i](X);

        return u;
    }

    std::array<autodiff::Vector3dual, NumBasisFunctions> uis;
    std::array<basis_function_op_type, NumBasisFunctions> phis;
};

template <class BasisFunctionType>
struct interpolation_op_t<BasisFunctionType, -1>
{
    using basis_function_op_type = BasisFunctionType;

    interpolation_op_t(
        std::vector<autodiff::Vector3dual> const& uis,
        std::vector<basis_function_op_type> const& phis)
        : uis(uis), phis(phis)
    {
    }

    autodiff::Vector3dual operator()(autodiff::Vector3dual const& X) const
    {
        autodiff::Vector3dual u;
        u << 0., 0., 0.;

        for (auto i = 0u; i < uis.size(); ++i)
            u += uis[i] * phis[i](X);

        return u;
    }

    std::vector<autodiff::Vector3dual> uis;
    std::vector<basis_function_op_type> phis;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_INTERPOLATION_H