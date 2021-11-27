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

namespace differentiable {

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
            kernel_type const& Wj                              = Wjs[k];
            Eigen::Vector<autodiff::dual, num_coefficients> Pj = polynomial3d<Order>(X - Xj);
            autodiff::dual WjX                                 = Wj(X);
            Eigen::Matrix<autodiff::dual, num_coefficients, num_coefficients> PjPjT =
                Pj * Pj.transpose();
            moment_matrix_type Ainc = WjX * Pj * Pj.transpose();
            // moment matrix gradient will be singular when Xj == Xi
            if (Xj.isApprox(X))
            {
                for (auto i = 0; i < num_coefficients; ++i)
                {
                    for (auto j = 0; j < num_coefficients; ++j)
                    {
                        Ainc(i, j).grad = 0.;
                    }
                }
            }
            A += Ainc; // W(X-Xj) P(Xj) P(Xj)^T
        }

        Eigen::Vector<autodiff::dual, num_coefficients> P  = polynomial3d<Order>(X - X);
        Eigen::Vector<autodiff::dual, num_coefficients> Pi = polynomial3d<Order>(X - Xi);

#if defined(_DEBUG)
        scalar_type const det = A.cast<scalar_type>().determinant();
        bool const invertible = std::abs(det) > sbs::eps();
        assert(invertible);
#endif
        moment_matrix_type Ainv = A.inverse();
        autodiff::dual WiX      = Wi(X);
        if (X.isApprox(Xi))
        {
            WiX.grad = 0.;
        }
        Eigen::Vector<autodiff::dual, num_coefficients> WiP     = WiX * P;
        Eigen::Vector<autodiff::dual, num_coefficients> AinvWiP = Ainv * WiP;
        // phi = P^T A(X)^-1 W(X-Xi) P(Xi)
        autodiff::dual phi = Pi.dot(AinvWiP);
        return phi;
    }

    autodiff::Vector3dual Xi;
    kernel_type Wi;
    std::vector<autodiff::Vector3dual> Xjs;
    std::vector<kernel_type> Wjs;
};

template <class KernelType = poly6_kernel_t>
struct sph_basis_function_t
{
    using kernel_type      = KernelType;
    sph_basis_function_t() = default;
    sph_basis_function_t(kernel_type const& Wi, scalar_type const Vi) : Wi(Wi), Vi(Vi) {}

    autodiff::dual operator()(autodiff::Vector3dual X) const
    {
        auto const Wij = Wi(X);
        return Vi * Wij;
    }

    kernel_type Wi;
    autodiff::dual Vi;
};

} // namespace differentiable

template <unsigned int Order>
struct polynomial_hat_basis_function_t
{
    using coefficients_type =
        Eigen::Vector<scalar_type, num_nodes_for_tetrahedron_cell_of_order_t<Order>::value>;
    static_assert(
        num_nodes_for_tetrahedron_cell_of_order_t<Order>::value ==
            num_coefficients_for_polynomial_of_order_t<Order>::value,
        "Number of nodes required in the tetrahedron should be same number of coefficients "
        "required in polynomial interpolation");

    polynomial_hat_basis_function_t() = default;
    polynomial_hat_basis_function_t(coefficients_type const& a) : a(a) {}

    scalar_type operator()(Eigen::Vector3d const& X) const { return eval(X); }

    scalar_type eval(Eigen::Vector3d const& X) const
    {
        auto const p = polynomial3d<Order>(X);
        return p.dot(a);
    }

    Eigen::Vector3d grad(Eigen::Vector3d const& X) const
    {
        static_assert(
            Order == 1u,
            "Gradient of higher-order polynomial functions not yet implemented");

        // a0 + a1*X + a2*Y + a3*Z;
        return Eigen::Vector3d{a(1), a(2), a(3)};
    }

    coefficients_type a;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_BASIS_FUNCTIONS_H
