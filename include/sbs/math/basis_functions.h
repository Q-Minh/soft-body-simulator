#ifndef SBS_MATH_BASIS_FUNCTIONS_H
#define SBS_MATH_BASIS_FUNCTIONS_H

#include "common.h"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace sbs {
namespace math {

template <unsigned int Order>
constexpr unsigned int num_nodes_for_tetrahedron_cell_of_order()
{
    unsigned int num_nodes                     = 0u;
    unsigned int constexpr num_samples_on_edge = Order + 1u;
    for (auto k = 0u; k < num_samples_on_edge; ++k)
    {
        auto const j_samples = num_samples_on_edge - k;
        for (auto j = 0u; j < j_samples; ++j)
        {
            auto const i_samples = num_samples_on_edge - k - j;
            for (auto i = 0u; i < i_samples; ++i)
            {
                ++num_nodes;
            }
        }
    }
    return num_nodes;
}

template <unsigned int Order>
struct num_nodes_for_tetrahedron_cell_of_order_t
{
    static constexpr unsigned int value = num_nodes_for_tetrahedron_cell_of_order<Order>();
};

struct linear_hat_basis_function_op_t
{
    linear_hat_basis_function_op_t(autodiff::Vector4dual const& a) : a(a) {}

    autodiff::dual operator()(autodiff::Vector3dual const& X) const
    {
        autodiff::Vector4dual p1 = polynomial3d<1>(X);
        return p1.dot(a);
    }

    autodiff::Vector4dual a;
};

template <unsigned int Order>
struct polynomial_hat_basis_function_t
{
    using coefficients_type =
        Eigen::Vector<autodiff::dual, num_nodes_for_tetrahedron_cell_of_order_t<Order>::value>;
    static_assert(
        num_nodes_for_tetrahedron_cell_of_order_t<Order>::value ==
            num_coefficients_for_polynomial_of_order_t > Order > ::value_compare,
        "Number of nodes required in the tetrahedron should be same number of coefficients "
        "required in polynomial interpolation");

    polynomial_hat_basis_function_t(coefficients_type const& a) : a(a) {}

    autodiff::dual operator()(autodiff::Vector3dual const& X) const
    {
        auto const p1 = polynomial3d<Order>(X);
        return p1.dot(a);
    }

    coefficients_type a;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_BASIS_FUNCTIONS_H
