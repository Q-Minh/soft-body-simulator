#include <sbs/physics/body.h>
#include <sbs/physics/collision/cd_system.h>
#include <sbs/physics/collision/collision_model.h>
#include <sbs/physics/simulation.h>
#include <sbs/physics/solver.h>
#include <sbs/physics/timestep.h>

namespace sbs {
namespace physics {

void timestep_t::step(simulation_t& simulation)
{
    scalar_type const dt = dt_ / static_cast<scalar_type>(substeps_);

    auto& particles = simulation.particles();

    // move particles using semi-implicit integration
    for (std::vector<particle_t>& body_particles : particles)
    {
        for (particle_t& p : body_particles)
        {
            p.v()  = p.v() + p.a() * dt;
            p.xi() = p.x() + p.v() * dt;
        }
    }

    // TODO: Cut
    // cut(simulation);
    // body->update_physical_model();

    auto const& cd_system = simulation.collision_detection_system();
    cd_system->execute();

    for (std::size_t s = 0u; s < substeps_; ++s)
    {
        solver_->solve(simulation, dt, iterations_);
    }

    simulation.collision_constraints().clear();

    auto& bodies = simulation.bodies();
    for (auto& body : bodies)
    {
        body->update_visual_model(simulation);
        body->update_collision_model(simulation);
    }

    cd_system->update(simulation);
}

} // namespace physics
} // namespace sbs