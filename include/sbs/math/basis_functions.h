#ifndef SBS_MATH_BASIS_FUNCTIONS_H
#define SBS_MATH_BASIS_FUNCTIONS_H

#include "common.h"
#include "kernels.h"

#include <Eigen/Geometry>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <vector>

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

    autodiff::dual operator()(autodiff::Vector3dual X) const
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
            num_coefficients_for_polynomial_of_order_t<Order>::value,
        "Number of nodes required in the tetrahedron should be same number of coefficients "
        "required in polynomial interpolation");

    polynomial_hat_basis_function_t() = default;
    polynomial_hat_basis_function_t(coefficients_type const& a) : a(a) {}

    autodiff::dual operator()(autodiff::Vector3dual X) const
    {
        auto const p1 = polynomial3d<Order>(X);
        return p1.dot(a);
    }

    coefficients_type a;
};

template <class KernelType = quartic_spline_kernel_t, unsigned int Order = 1>
struct mls_basis_function_t
{
    using kernel_type        = KernelType;
    using moment_matrix_type = Eigen::Matrix<
        autodiff::dual,
        num_coefficients_for_polynomial_of_order_t<Order>::value,
        num_coefficients_for_polynomial_of_order_t<Order>::value>;

    mls_basis_function_t() = default;
    mls_basis_function_t(
        autodiff::Vector3dual Xi,
        kernel_type Wi,
        std::vector<autodiff::Vector3dual> const& Xjs,
        std::vector<kernel_type> const& Wjs)
        : Xi(Xi), Wi(Wi), Xjs(Xjs), Wjs(Wjs)
    {
        assert(Xjs.size() == Wjs.size());
    }

    autodiff::dual operator()(autodiff::Vector3dual X) const
    {
        // moment matrix
        moment_matrix_type A{};
        A.setZero();

        unsigned int constexpr num_coefficients =
            num_coefficients_for_polynomial_of_order_t<Order>::value;

        for (auto k = 0u; k < Xjs.size(); ++k)
        {
            autodiff::Vector3dual const& Xj                    = Xjs[k];
            quartic_spline_kernel_t const& Wj                  = Wjs[k];
            Eigen::Vector<autodiff::dual, num_coefficients> P1 = polynomial3d<Order>(Xj);
            A += Wj(X) * P1 * P1.transpose(); // W(X-Xj) P(Xj) P(Xj)^T
        }

        Eigen::Vector<autodiff::dual, num_coefficients> P        = polynomial3d<Order>(X);
        Eigen::Vector<autodiff::dual, num_coefficients> Pi       = polynomial3d<Order>(Xi);
        moment_matrix_type Ainv                                  = A.inverse();
        Eigen::Vector<autodiff::dual, num_coefficients> WiPi     = Wi(X) * Pi;
        Eigen::Vector<autodiff::dual, num_coefficients> AinvWiPi = Ainv * WiPi;
        // phi = P^T A(X)^-1 W(X-Xi) P(Xi)
        autodiff::dual phi = P.dot(AinvWiPi);
        return phi;
    }

    autodiff::Vector3dual Xi;
    kernel_type Wi;
    std::vector<autodiff::Vector3dual> Xjs;
    std::vector<kernel_type> Wjs;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_BASIS_FUNCTIONS_H
