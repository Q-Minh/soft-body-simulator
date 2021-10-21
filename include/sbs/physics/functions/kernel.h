#ifndef SBS_PHYSICS_FUNCTIONS_KERNEL_H
#define SBS_PHYSICS_FUNCTIONS_KERNEL_H

#include <Eigen/Core>
#include <sbs/aliases.h>

namespace sbs {
namespace physics {
namespace functions {

class kernel_t
{
  public:
    kernel_t() = default;
    kernel_t(Eigen::Vector3d const& xi, scalar_type const h);

    scalar_type h() const;
    scalar_type& h();
    Eigen::Vector3d const& xi() const;
    Eigen::Vector3d& xi();

    virtual scalar_type operator()(Eigen::Vector3d const& xj) const = 0;
    virtual Eigen::Vector3d grad(Eigen::Vector3d const& xj) const   = 0;

  private:
    Eigen::Vector3d xi_;
    scalar_type h_;
};

class poly6_kernel_t : public kernel_t
{
  public:
    poly6_kernel_t() = default;
    poly6_kernel_t(Eigen::Vector3d const& xi, scalar_type const h);

    poly6_kernel_t(poly6_kernel_t const& other) = default;
    poly6_kernel_t(poly6_kernel_t&& other)      = default;

    poly6_kernel_t& operator=(poly6_kernel_t const& other) = default;
    poly6_kernel_t& operator=(poly6_kernel_t&& other) = default;

    virtual scalar_type operator()(Eigen::Vector3d const& xj) const override;
    virtual Eigen::Vector3d grad(Eigen::Vector3d const& xj) const override;

  private:
    static scalar_type constexpr mult_315_64_pi = 315. / (64. * 3.14159265359);

    scalar_type h2_;
    scalar_type h9_;
};

} // namespace functions
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_FUNCTIONS_KERNEL_H