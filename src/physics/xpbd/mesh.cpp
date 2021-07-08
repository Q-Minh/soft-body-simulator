#include "physics/xpbd/mesh.h"

namespace sbs {
namespace physics {
namespace xpbd {

tetrahedral_mesh_t::tetrahedral_mesh_t(
    common::geometry_t const& geometry,
    simulation_parameters_t const& simulation_params,
    build_topology_parameters_t const& topology_params)
    : physics::renderable_topological_simulated_tetrahedral_mesh_t{geometry, topology_params},
      simulation_params_{simulation_params}
{
}

simulation_parameters_t const& sbs::physics::xpbd::tetrahedral_mesh_t::simulation_parameters() const
{
    return simulation_params_;
}

simulation_parameters_t& sbs::physics::xpbd::tetrahedral_mesh_t::simulation_parameters()
{
    return simulation_params_;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs