#ifndef SBS_PHYSICS_XPBD_MESHLESS_COROTATIONAL_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_MESHLESS_COROTATIONAL_CONSTRAINT_H

#include <sbs/aliases.h>
#include <sbs/physics/constraint.h>
#include <sbs/physics/mechanics/meshless_node.h>

namespace sbs {
namespace physics {
namespace xpbd {

class meshless_stvk_constraint_t : public constraint_t
{
  public:
    meshless_stvk_constraint_t(
        scalar_type const alpha,
        scalar_type const beta,
        simulation_t const& simulation,
        index_type bi,
        index_type vi,
        scalar_type young_modulus,
        scalar_type poisson_ratio,
        mechanics::meshless_node_t const& node);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    mechanics::meshless_node_t const& node_;
    index_type bi_;
    index_type vi_;
    scalar_type mu_;
    scalar_type lambda_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_MESHLESS_COROTATIONAL_CONSTRAINT_H