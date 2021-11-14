#ifndef SBS_PHYSICS_XPBD_PARTICLE_H
#define SBS_PHYSICS_XPBD_PARTICLE_H

#include "sbs/aliases.h"

#include <Eigen/Core>

namespace sbs {
namespace physics {
namespace xpbd {

class particle_t
{
  public:
    using position_type     = Eigen::Vector3d;
    using velocity_type     = Eigen::Vector3d;
    using acceleration_type = Eigen::Vector3d;
    using force_type        = Eigen::Vector3d;

    particle_t() = default;
    particle_t(position_type const& p);

    particle_t(particle_t const& other) = default;
    particle_t(particle_t&& other)      = default;
    particle_t& operator=(particle_t const& other) = default;
    particle_t& operator=(particle_t&& other) = default;

    position_type const& x0() const;
    position_type const& x() const;
    position_type const& xi() const;
    position_type const& xn() const;
    velocity_type const& v() const;
    force_type const& f() const;
    scalar_type const& mass() const;
    scalar_type invmass() const;
    bool fixed() const;
    acceleration_type a() const;

    position_type& x0();
    position_type& x();
    position_type& xi();
    position_type& xn();
    velocity_type& v();
    force_type& f();
    scalar_type& mass();

  private:
    position_type x0_;
    position_type xi_;
    position_type xn_;
    position_type x_;
    velocity_type v_;
    force_type f_;
    scalar_type m_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_PARTICLE_H