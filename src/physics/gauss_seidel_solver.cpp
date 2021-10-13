#include <sbs/physics/constraint.h>
#include <sbs/physics/gauss_seidel_solver.h>
#include <sbs/physics/simulation.h>

namespace sbs {
namespace physics {

void gauss_seidel_solver_t::solve(simulation_t& simulation, scalar_type dt, std::size_t iterations)
{
    auto& particles = simulation.particles();

    auto& collision_constraints = simulation.collision_constraints();
    auto& constraints           = simulation.constraints();

    for (auto& collision_constraint : collision_constraints)
    {
        collision_constraint->prepare_for_projection(simulation);
    }
    for (auto& constraint : constraints)
    {
        constraint->prepare_for_projection(simulation);
    }

    // solver loop
    for (std::size_t k = 0u; k < iterations; ++k)
    {
        // solve positional constraints
        for (auto& collision_constraint : collision_constraints)
        {
            collision_constraint->project_positions(simulation, dt);
        }
        for (auto& constraint : constraints)
        {
            constraint->project_positions(simulation, dt);
        }
    }
}

} // namespace physics
} // namespace sbs