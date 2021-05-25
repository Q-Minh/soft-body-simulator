#ifndef SBS_PHYSICS_XPBD_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_CONSTRAINT_H

#include "common/scene.h"

#include <cstdint>
#include <utility>

namespace sbs {
namespace physics {
namespace xpbd {

class constraint_t
{
  public:
    using scalar_type = double;
    using index_type  = std::pair<std::uint32_t, std::uint32_t>;

    constraint_t(scalar_type const alpha) : alpha_(alpha) {}

    virtual void project(
        common::scene_t& scene,
        scalar_type& lagrange_multiplier,
        scalar_type const dt) const = 0;

  protected:
    scalar_type alpha_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_CONSTRAINT_H