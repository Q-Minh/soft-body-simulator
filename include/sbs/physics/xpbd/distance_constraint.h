#ifndef SBS_PHYSICS_XPBD_DISTANCE_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_DISTANCE_CONSTRAINT_H

#include "sbs/aliases.h"
#include "sbs/physics/xpbd/constraint.h"

#include <Eigen/Core>

namespace sbs {
namespace physics {
namespace xpbd {

class distance_constraint_t : public constraint_t
{
  public:
    distance_constraint_t(
        scalar_type const alpha,
        scalar_type const beta,
        physics::xpbd::simulation_t const& simulation,
        index_type b1,
        index_type b2,
        index_type v1,
        index_type v2);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    index_type b1_;
    index_type v1_;
    index_type b2_;
    index_type v2_;
    scalar_type d_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_DISTANCE_CONSTRAINT_H