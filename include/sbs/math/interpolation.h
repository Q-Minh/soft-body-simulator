#ifndef SBS_MATH_INTERPOLATION_H
#define SBS_MATH_INTERPOLATION_H

#include <array>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <vector>

namespace sbs {
namespace math {

enum interpolation_op_variant { dynamic_owning = -1, dynamic_non_owning = -2 };

namespace differentiable {

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

template <class BasisFunctionType>
struct interpolation_op_t<BasisFunctionType, -2>
{
    using basis_function_op_type = BasisFunctionType;

    interpolation_op_t(
        std::vector<autodiff::Vector3dual> const& uis,
        std::vector<basis_function_op_type> const& phis,
        std::vector<index_type> const& is)
        : uis(uis), phis(phis), is(is)
    {
    }

    autodiff::Vector3dual operator()(autodiff::Vector3dual const& X) const
    {
        autodiff::Vector3dual u;
        u << 0., 0., 0.;

        for (auto i = 0u; i < is.size(); ++i)
            u += uis[is[i]] * phis[is[i]](X);

        return u;
    }

    std::vector<autodiff::Vector3dual> const* uis;
    std::vector<basis_function_op_type> const* phis;
    std::vector<index_type> is;
};

} // namespace differentiable
} // namespace math
} // namespace sbs

#endif // SBS_MATH_INTERPOLATION_H