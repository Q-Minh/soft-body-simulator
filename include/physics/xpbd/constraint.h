#ifndef SBS_PHYSICS_XPBD_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_CONSTRAINT_H

#include "common/scene.h"

#include <Eigen/Core>
#include <cstdint>

namespace sbs {
namespace physics {
namespace xpbd {

class constraint_t
{
  public:
    using scalar_type    = double;
    using index_type     = std::uint32_t;
    using positions_type = Eigen::Matrix3Xd;
    using masses_type    = Eigen::VectorXd;

    constraint_t(scalar_type const alpha) : alpha_(alpha) {}

    virtual void project(
        std::vector<positions_type>& positions,
        std::vector<masses_type> const& masses,
        scalar_type& lagrange_multiplier,
        scalar_type const dt) const = 0;

  protected:
    scalar_type alpha_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_CONSTRAINT_H