#ifndef SBS_PHYSICS_XPBD_COLLISION_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_COLLISION_CONSTRAINT_H

#include "sbs/physics/xpbd/constraint.h"

#include <Eigen/Core>

namespace sbs {
namespace physics {
namespace xpbd {

// Forward declares
class simulation_t;

class collision_constraint_t : public constraint_t
{
  public:
    collision_constraint_t(
        scalar_type alpha /*compliance*/,
        scalar_type beta /*damping*/,
        simulation_t const& simulation,
        index_type bi /*penetrating body*/,
        index_type vi /*penetrating vertex*/,
        Eigen::Vector3d const& p, /*contact point*/
        Eigen::Vector3d const& n /*surface normal of correction*/);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;
    scalar_type evaluate(Eigen::Vector3d const& p) const;

  private:
    index_type bi_; ///< Penetrating body
    index_type vi_; ///< Index of penetrating vertex

    Eigen::Vector3d qs_; ///< Intersection point
    Eigen::Vector3d n_;  ///< Normal at intersection point
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_COLLISION_CONSTRAINT_H