#include <sbs/physics/functions/kernel.h>

namespace sbs {
namespace physics {
namespace functions {

kernel_t::kernel_t(Eigen::Vector3d const& xi, scalar_type const h) : xi_(xi), h_(h) {}

scalar_type kernel_t::h() const
{
    return h_;
}

Eigen::Vector3d const& kernel_t::xi() const
{
    return xi_;
}

poly6_kernel_t::poly6_kernel_t(Eigen::Vector3d const& xi, scalar_type const h)
    : kernel_t(xi, h), h2_(h * h), h9_(h * h * h * h * h * h * h * h * h)
{
}

scalar_type poly6_kernel_t::operator()(Eigen::Vector3d const& xj) const
{
    Eigen::Vector3d const diff = xj - xi();
    scalar_type const r2       = diff.squaredNorm();

    if (r2 > h2_)
        return 0.;

    auto const mult       = mult_315_64_pi / h9_;
    scalar_type const hr3 = (h2_ - r2) * (h2_ - r2) * (h2_ - r2);
    return mult * hr3;
}

Eigen::Vector3d poly6_kernel_t::grad(Eigen::Vector3d const& xj) const
{
    Eigen::Vector3d const diff  = xj - xi();
    scalar_type const r         = diff.norm();
    scalar_type const r2        = r * r;
    scalar_type const mult      = -6. * r * mult_315_64_pi / h9_;
    scalar_type const hr3       = (h2_ - r2) * (h2_ - r2) * (h2_ - r2);
    scalar_type const dWdr      = mult * hr3;
    Eigen::Vector3d const drdxj = diff.normalized();
    return dWdr * drdxj;
}

} // namespace functions
} // namespace physics
} // namespace sbs