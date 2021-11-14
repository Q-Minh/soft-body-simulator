#include "sbs/physics/timestep.h"

#include "sbs/physics/body/body.h"
#include "sbs/physics/collision/cd_system.h"
#include "sbs/physics/collision/collision_model.h"
#include "sbs/physics/xpbd/simulation.h"
#include "sbs/physics/xpbd/solver.h"

#include <execution>

namespace sbs {
namespace physics {

timestep_t::timestep_t(
    scalar_type const dt,
    std::size_t const iterations,
    std::size_t const substeps)
    : dt_(dt), iterations_(iterations), substeps_(substeps)
{
}

void timestep_t::step(xpbd::simulation_t& simulation)
{
    scalar_type const dt = dt_ / static_cast<scalar_type>(substeps_);

    auto& particles = simulation.particles();

    // TODO: Cut
    // cut(simulation);
    // body->update_physical_model();
    auto const& cd_system = simulation.collision_detection_system();
    cd_system->contact_handler()->on_cd_starting();
    cd_system->execute();
    cd_system->contact_handler()->on_cd_ending();

    // if (!simulation.collision_constraints().empty())
    //{
    //     std::cout << "Collisions: " << simulation.collision_constraints().size() << "\n";
    // }

    for (std::size_t s = 0u; s < substeps_; ++s)
    {
        // move particles using semi-implicit integration
        for (std::vector<xpbd::particle_t>& body_particles : particles)
        {
            std::for_each(body_particles.begin(), body_particles.end(), [dt](xpbd::particle_t& p) {
                p.f().y() -= scalar_type{9.81};
                p.v()  = p.v() + p.a() * dt;
                p.xi() = p.x() + p.v() * dt;
            });
        }

        solver_->solve(simulation, dt, iterations_);

        // set solution
        for (std::vector<xpbd::particle_t>& body_particles : particles)
        {
            std::for_each(body_particles.begin(), body_particles.end(), [dt](xpbd::particle_t& p) {
                p.x()  = p.xi();
                p.v()  = (p.x() - p.xn()) / dt;
                p.xn() = p.x();
                p.f().setZero();
            });
        }
    }

    auto& bodies = simulation.bodies();
    for (auto& body : bodies)
    {
        body->update_physical_model();
        body->update_visual_model();
        body->update_collision_model();
        body->visual_model().mark_vertices_dirty();
    }

    simulation.collision_constraints().clear();
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

std::unique_ptr<xpbd::solver_t> const& timestep_t::solver() const
{
    return solver_;
}
std::unique_ptr<xpbd::solver_t>& timestep_t::solver()
{
    return solver_;
}

} // namespace physics
} // namespace sbs