#ifndef SBS_PHYSICS_XPBD_COLLISION_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_COLLISION_CONSTRAINT_H

#include "constraint.h"

namespace sbs {
namespace physics {
namespace xpbd {

class collision_constraint_t : public constraint_t
{
  public:
    using scalar_type     = typename constraint_t::scalar_type;
    using index_type      = typename constraint_t::index_type;
    using index_pair_type = std::pair<index_type, index_type>;
    using positions_type  = typename constraint_t::positions_type;
    using position_type   = Eigen::Vector3d;
    using normal_type     = Eigen::Vector3d;
    using masses_type     = typename constraint_t::masses_type;

    collision_constraint_t(
        //scalar_type const alpha,
        index_pair_type penetrating_vertex,
        position_type const& q,
        normal_type const& n);

    virtual void project(
        std::vector<positions_type>& positions,
        std::vector<masses_type> const& masses,
        scalar_type& lagrange_multiplier,
        scalar_type const dt) const override;

    scalar_type evaluate(position_type const& p) const;

  private:
    index_type p_;  ///< Penetrating vertex
    index_type b1_; ///< The body that owns vertex q

    position_type q_; ///< Intersection point
    normal_type n_;   ///< Normal at intersection point
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_COLLISION_CONSTRAINT_H