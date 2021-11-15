#ifndef SBS_MATH_KERNELS_H
#define SBS_MATH_KERNELS_H

#include "sbs/aliases.h"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace sbs {
namespace math {

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
        autodiff::dual r = (X - Xi).norm();
        if (r > h)
            return 0.;

        autodiff::dual r2   = r * r;
        autodiff::dual h_r  = (h2 - r2);
        autodiff::dual h_r3 = h_r * h_r * h_r;
        return alpha * h_r3;
    }

    static constexpr scalar_type pi = 3.1415926535;

    autodiff::Vector3dual Xi;
    autodiff::dual h;
    autodiff::dual h2;
    autodiff::dual h9;
    autodiff::dual alpha;
};

struct quartic_spline_kernel_t
{
    quartic_spline_kernel_t() = default;
    quartic_spline_kernel_t(autodiff::Vector3dual Xi, scalar_type r) : Xi(Xi), r(r) {}

    autodiff::dual operator()(autodiff::Vector3dual X) const
    {
        autodiff::dual d  = (X - Xi).norm() / r;
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