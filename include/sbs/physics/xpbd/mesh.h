#ifndef SBS_PHYSICS_XPBD_MESH_H
#define SBS_PHYSICS_XPBD_MESH_H

#include "sbs/physics/topological_mesh.h"
#include "sbs/physics/xpbd/simulation_parameters.h"

namespace sbs {
namespace physics {
namespace xpbd {

class tetrahedral_mesh_t : public physics::renderable_topological_simulated_tetrahedral_mesh_t
{
  public:
    tetrahedral_mesh_t(
        common::geometry_t const& geometry,
        simulation_parameters_t const& simulation_params,
        build_topology_parameters_t const& topology_params = build_topology_parameters_t{});

    simulation_parameters_t const& simulation_parameters() const;
    simulation_parameters_t& simulation_parameters();

  private:
    simulation_parameters_t simulation_params_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_MESH_H