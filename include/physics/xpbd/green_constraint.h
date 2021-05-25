#ifndef SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H

#include "constraint.h"

namespace sbs {
namespace physics {
namespace xpbd {

class green_constraint_t : public constraint_t
{
  public:
    using scalar_type = typename constraint_t::scalar_type;
    using index_type  = typename constraint_t::index_type;

    green_constraint_t(
        scalar_type const alpha,
        index_type const& v1,
        index_type const& v2,
        index_type const& v3,
        index_type const& v4)
        : constraint_t(alpha), v1_(v1), v2_(v2), v3_(v3), v4_(v4)
    {
    }

    virtual void project(
        common::scene_t& scene,
        scalar_type& lagrange_multiplier,
        scalar_type const dt) const override;

  private:
    index_type v1_;
    index_type v2_;
    index_type v3_;
    index_type v4_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H