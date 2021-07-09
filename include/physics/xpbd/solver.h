#ifndef SBS_PHYSICS_XPBD_SOLVER_H
#define SBS_PHYSICS_XPBD_SOLVER_H

#include "constraint.h"

#include <Eigen/Core>
#include <array>
#include <map>
#include <memory>
#include <vector>

namespace sbs {

// forward declares
namespace common {

class renderable_node_t;
class shared_vertex_surface_mesh_i;

} // namespace common

namespace physics {
namespace xpbd {

// forward declare
class tetrahedral_mesh_t;

class solver_t
{
  public:
    solver_t() = default;
    solver_t(double timestep, std::uint32_t iterations, std::uint32_t substeps);
    solver_t(std::vector<std::shared_ptr<common::renderable_node_t>> const& bodies);
    solver_t(
        double timestep,
        std::uint32_t iterations,
        std::uint32_t substeps,
        std::vector<std::shared_ptr<common::renderable_node_t>> const& bodies);

    void setup(std::vector<std::shared_ptr<common::renderable_node_t>> const& bodies);
    void reset();
    void step();

    double const& timestep() const;
    double& timestep();

    std::uint32_t const& iterations() const;
    std::uint32_t& iterations();

    std::uint32_t const& substeps() const;
    std::uint32_t& substeps();

    std::vector<std::shared_ptr<xpbd::tetrahedral_mesh_t>> const& simulated_bodies() const;

  protected:
    void handle_collisions(std::vector<Eigen::Matrix3Xd> const& P);

    void create_green_constraints_for_body(xpbd::tetrahedral_mesh_t* body);
    void create_distance_constraints_for_body(xpbd::tetrahedral_mesh_t* body);

  private:
    using constraint_map_key_type = std::pair<xpbd::tetrahedral_mesh_t*, std::uint32_t>;

    struct constraint_map_key_less
    {
        bool
        operator()(constraint_map_key_type const& key1, constraint_map_key_type const& key2) const
        {
            std::pair<std::uintptr_t, std::uint32_t> const a1{
                reinterpret_cast<std::uintptr_t>(key1.first),
                key1.second};

            std::pair<std::uintptr_t, std::uint32_t> const a2{
                reinterpret_cast<std::uintptr_t>(key2.first),
                key2.second};

            return a1 < a2;
        }
    };

    std::vector<std::shared_ptr<xpbd::tetrahedral_mesh_t>>
        physics_bodies_; ///< List of XPBD based physically simulated bodies
    std::vector<std::shared_ptr<common::shared_vertex_surface_mesh_i>>
        environment_bodies_; ///< List of non-simulated bodies participating in collision
                             ///< detection/handling

    std::vector<std::unique_ptr<constraint_t>> constraints_; ///< List of XPBD constraints
    std::vector<std::unique_ptr<constraint_t>>
        collision_constraints_; ///< List of XPBD collision constraints that is rebuilt every frame

    std::map<
        constraint_map_key_type,
        std::uint32_t /* index of constraint */,
        constraint_map_key_less>
        tetrahedron_to_constraint_map_; ///< Mapping from body/tetrahedron to its corresponding
                                        ///< constraint's index in the list of constraints

    double dt_;                                    ///< Timestep duration of the simulation
    std::uint32_t substeps_;                       ///< Number of substeps in one solver step
    std::uint32_t iteration_count_;                ///< Number of iterations per substep
    std::vector<std::vector<Eigen::Vector3d>> x0_; ///< Initial positions of one timestep
    std::vector<double> lagrange_multipliers_;     ///< Lagrange multipliers of XPBD formulation
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_SOLVER_H