#include <Eigen/Geometry>
#include <sbs/physics/simulation.h>
#include <sbs/physics/xpbd/collision_constraint.h>

namespace sbs {
namespace physics {
namespace xpbd {

collision_constraint_t::collision_constraint_t(
    scalar_type alpha,
    scalar_type beta,
    simulation_t const& simulation,
    index_type bi,
    index_type vi,
    Eigen::Vector3d const& p,
    Eigen::Vector3d const& n)
    : constraint_t{alpha, beta}, bi_(bi), vi_(vi), qs_(p), n_(n)
{
}

void collision_constraint_t::project_positions(simulation_t& simulation, scalar_type dt)
{
    particle_t& p       = simulation.particles()[bi_][vi_];
    scalar_type const w = p.invmass();
    scalar_type const C = evaluate(p.xi());

    if (C >= static_cast<scalar_type>(0.))
        return;

    scalar_type const alpha_tilde    = alpha_ / (dt * dt);
    scalar_type const delta_lagrange = -(C + alpha_tilde * lagrange_) / (w + alpha_tilde);

    lagrange_ += delta_lagrange;

    /**
     * grad(C) = n
     * ||n||^2 = 1,
     * so w*||n||^2 = w
     */
    p.xi() += w * n_ * delta_lagrange;
}

scalar_type collision_constraint_t::evaluate(Eigen::Vector3d const& p) const
{
    Eigen::Vector3d const qp = p - qs_;
    double const C           = qp.dot(n_);
    return C;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs
