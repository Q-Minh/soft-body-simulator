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

    // TODO: Cut
    // cut(simulation);
    // body->update_physical_model();

    auto const& cd_system = simulation.collision_detection_system();
    cd_system->execute();

    for (std::size_t s = 0u; s < substeps_; ++s)
    {
        // move particles using semi-implicit integration
        for (std::vector<particle_t>& body_particles : particles)
        {
            for (particle_t& p : body_particles)
            {
                p.f().y() -= scalar_type{9.81};
                p.v()  = p.v() + p.a() * dt;
                p.xi() = p.x() + p.v() * dt;
            }
        }

        solver_->solve(simulation, dt, iterations_);

        // set solution
        for (std::vector<particle_t>& body_particles : particles)
        {
            for (particle_t& p : body_particles)
            {
                p.x()  = p.xi();
                p.v()  = (p.x() - p.xn()) / dt;
                p.xn() = p.x();
                p.f().setZero();
            }
        }
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

scalar_type timestep_t::dt() const
{
    return dt_;
}
scalar_type& timestep_t::dt()
{
    return dt_;
}

std::size_t timestep_t::iterations() const
{
    return iterations_;
}
std::size_t& timestep_t::iterations()
{
    return iterations_;
}

std::size_t timestep_t::substeps() const
{
    return substeps_;
}
std::size_t& timestep_t::substeps()
{
    return substeps_;
}

std::unique_ptr<solver_t> const& timestep_t::solver() const
{
    return solver_;
}
std::unique_ptr<solver_t>& timestep_t::solver()
{
    return solver_;
}

} // namespace physics
} // namespace sbs