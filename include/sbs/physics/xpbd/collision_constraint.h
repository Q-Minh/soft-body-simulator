#ifndef SBS_PHYSICS_XPBD_COLLISION_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_COLLISION_CONSTRAINT_H

#include "sbs/physics/xpbd/constraint.h"

#include <Eigen/Core>

namespace sbs {

// forward declares
namespace common {

struct triangle_t;

} // namespace common

namespace physics {
namespace xpbd {

class collision_constraint_t : public constraint_t
{
  public:
    using scalar_type   = double;
    using position_type = Eigen::Vector3d;
    using normal_type   = Eigen::Vector3d;
    using index_type    = std::uint32_t;
    using body_ptr_type = tetrahedral_mesh_t*;

    collision_constraint_t(
        scalar_type alpha /*compliance*/,
        body_ptr_type b /*penetrating body*/,
        index_type vi /*penetrating vertex*/,
        common::triangle_t const& t /*penetrated triangle*/,
        normal_type const& n /*surface normal at qs*/);

    void project(
        std::vector<std::shared_ptr<tetrahedral_mesh_t>> const& bodies,
        scalar_type& lagrange_multiplier,
        scalar_type const dt) const;

    scalar_type evaluate(position_type const& p) const;

  private:
    body_ptr_type b_; ///< Penetrating body
    index_type vi_;   ///< Index of penetrating vertex

    position_type qs_; ///< Intersection point
    normal_type n_;    ///< Normal at intersection point
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_COLLISION_CONSTRAINT_H