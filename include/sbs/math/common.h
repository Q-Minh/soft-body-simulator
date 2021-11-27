#ifndef SBS_MATH_COMMON_H
#define SBS_MATH_COMMON_H

#include "sbs/aliases.h"

#include <Eigen/Core>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace sbs {
namespace math {

template <unsigned int Order>
constexpr unsigned int num_coefficients_for_polynomial_of_order()
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
struct num_coefficients_for_polynomial_of_order_t
{
    static constexpr unsigned int value = num_coefficients_for_polynomial_of_order<Order>();
};

namespace differentiable {

template <unsigned int Order>
Eigen::Vector<autodiff::dual, num_coefficients_for_polynomial_of_order_t<Order>::value>
polynomial3d(autodiff::Vector3dual const& X)
{
    static_assert(Order < 4 && Order >= 0, "Polynomials of order 0,1,2,3 supported");

    unsigned int constexpr num_coefficients =
        num_coefficients_for_polynomial_of_order_t<Order>::value;

    if constexpr (Order == 0)
    {
        return Eigen::Vector<autodiff::dual, num_coefficients>(1.);
    }
    if constexpr (Order == 1)
    {
        return Eigen::Vector<autodiff::dual, num_coefficients>(1., X.x(), X.y(), X.z());
    }
    if constexpr (Order == 2)
    {
        auto const& x = X.x();
        auto const& y = X.y();
        auto const& z = X.z();
        return Eigen::Vector<autodiff::dual, num_coefficients>(
            1.,
            x,
            y,
            z,
            x * x,
            x * y,
            x * z,
            y * y,
            y * z,
            z * z);
    }
    if constexpr (Order == 3)
    {
        auto const& x = X.x();
        auto const& y = X.y();
        auto const& z = X.z();
        auto const x2 = x * x;
        auto const y2 = y * y;
        auto const z2 = z * z;
        auto const x3 = x2 * x;
        auto const y3 = y2 * y;
        auto const z3 = z2 * z;
        return Eigen::Vector<autodiff::dual, num_coefficients>(
            1.,
            x,
            y,
            z,
            x2,
            x * y,
            x * z,
            y2,
            y * z,
            z2,
            x3,
            x2 * y,
            x2 * z,
            x * y2,
            x * y * z,
            x * z2,
            y3,
            y2 * z,
            y * z2,
            z3);
    }
}

} // namespace differentiable

template <unsigned int Order>
Eigen::Vector<scalar_type, num_coefficients_for_polynomial_of_order_t<Order>::value>
polynomial3d(Eigen::Vector3d const& X)
{
    static_assert(Order < 4 && Order >= 0, "Polynomials of order 0,1,2,3 supported");

    unsigned int constexpr num_coefficients =
        num_coefficients_for_polynomial_of_order_t<Order>::value;

    if constexpr (Order == 0)
    {
        return Eigen::Vector<scalar_type, num_coefficients>(1.);
    }
    if constexpr (Order == 1)
    {
        return Eigen::Vector<scalar_type, num_coefficients>(1., X.x(), X.y(), X.z());
    }
    if constexpr (Order == 2)
    {
        auto const& x = X.x();
        auto const& y = X.y();
        auto const& z = X.z();
        return Eigen::Vector<scalar_type, num_coefficients>(
            1.,
            x,
            y,
            z,
            x * x,
            x * y,
            x * z,
            y * y,
            y * z,
            z * z);
    }
    if constexpr (Order == 3)
    {
        auto const& x = X.x();
        auto const& y = X.y();
        auto const& z = X.z();
        auto const x2 = x * x;
        auto const y2 = y * y;
        auto const z2 = z * z;
        auto const x3 = x2 * x;
        auto const y3 = y2 * y;
        auto const z3 = z2 * z;
        return Eigen::Vector<scalar_type, num_coefficients>(
            1.,
            x,
            y,
            z,
            x2,
            x * y,
            x * z,
            y2,
            y * z,
            z2,
            x3,
            x2 * y,
            x2 * z,
            x * y2,
            x * y * z,
            x * z2,
            y3,
            y2 * z,
            y * z2,
            z3);
    }
}

} // namespace math
} // namespace sbs

#endif // SBS_MATH_COMMON_H