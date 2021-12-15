#include "sbs/physics/xpbd/simulation.h"

#include "sbs/physics/body/body.h"

namespace sbs {
namespace physics {
namespace xpbd {

void simulation_t::use_collision_detection_system(std::unique_ptr<collision::cd_system_t> cd_system)
{
    cd_system_ = std::move(cd_system);
}

void simulation_t::add_particle(particle_t const& p, index_type const body_idx)
{
    particles_[body_idx].push_back(p);
}

index_type simulation_t::add_body(std::unique_ptr<body::body_t> body)
{
    index_type const bi = static_cast<index_type>(bodies_.size());
    bodies_.push_back(std::move(body));
    particles_.push_back({});
    return bi;
}

index_type simulation_t::add_body()
{
    index_type const bi = static_cast<index_type>(bodies_.size());
    bodies_.push_back({});
    particles_.push_back({});
    return bi;
}

void simulation_t::add_constraint(std::unique_ptr<constraint_t> constraint)
{
    constraints_.push_back(std::move(constraint));
}

void simulation_t::remove_constraint(index_type const constraint_idx)
{
    std::size_t const last = constraints_.size() - 1u;
    std::swap(constraints_[constraint_idx], constraints_[last]);
    constraints_.erase(constraints_.begin() + last, constraints_.end());
}

void simulation_t::add_collision_constraint(std::unique_ptr<constraint_t> collision_constraint)
{
    collision_constraints_.push_back(std::move(collision_constraint));
}

void simulation_t::apply_dirichlet_boundary_conditions(
    std::vector<index_type> const& bs,
    std::vector<index_type> const& is,
    std::vector<Eigen::Vector3d> const& xis)
{
    assert(bs.size() == is.size());
    assert(bs.size() == xis.size());

    std::size_t const N = bs.size();
    for (std::size_t n = 0u; n < N; ++n)
    {
        index_type const b        = bs[n];
        index_type const i        = is[n];
        Eigen::Vector3d const& xi = xis[n];
        particle_t& pi            = particles_[b][i];
        pi.mass()                 = 0.;
        pi.x()                    = xi;
        pi.xn()                   = xi;
        pi.xi()                   = xi;
        pi.v().setZero();
    }
}

std::vector<std::vector<particle_t>> const& simulation_t::particles() const
{
    return particles_;
}
std::vector<std::vector<particle_t>>& simulation_t::particles()
{
    return particles_;
}
std::vector<std::unique_ptr<body::body_t>> const& simulation_t::bodies() const
{
    return bodies_;
}
std::vector<std::unique_ptr<body::body_t>>& simulation_t::bodies()
{
    return bodies_;
}
std::vector<std::unique_ptr<constraint_t>> const& simulation_t::constraints() const
{
    return constraints_;
}
std::vector<std::unique_ptr<constraint_t>>& simulation_t::constraints()
{
    return constraints_;
}
std::vector<std::unique_ptr<constraint_t>> const& simulation_t::collision_constraints() const
{
    return collision_constraints_;
}
std::vector<std::unique_ptr<constraint_t>>& simulation_t::collision_constraints()
{
    return collision_constraints_;
}
std::unique_ptr<collision::cd_system_t> const& simulation_t::collision_detection_system() const
{
    return cd_system_;
}
xpbd::simulation_parameters_t const& simulation_t::simulation_parameters() const
{
    return simulation_parameters_;
}
xpbd::simulation_parameters_t& simulation_t::simulation_parameters()
{
    return simulation_parameters_;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs