#include <sbs/physics/simulation.h>
#include <sbs/physics/xpbd/distance_constraint.h>

namespace sbs {
namespace physics {
namespace xpbd {

distance_constraint_t::distance_constraint_t(
    scalar_type const alpha,
    scalar_type const beta,
    simulation_t const& simulation,
    index_type b1,
    index_type b2,
    index_type v1,
    index_type v2)
    : constraint_t(alpha, beta), b1_(b1), b2_(b2), v1_(v1), v2_(v2)
{
    auto const& p1 = simulation.particles()[b1][v1];
    auto const& p2 = simulation.particles()[b2][v2];

    d_ = (p1.x0() - p2.x0()).norm();
}

void distance_constraint_t::project_positions(simulation_t& simulation, scalar_type dt)
{
    particle_t& p1 = simulation.particles()[b1_][v1_];
    particle_t& p2 = simulation.particles()[b2_][v2_];

    scalar_type const w1 = p1.invmass();
    scalar_type const w2 = p2.invmass();

    Eigen::Vector3d const diff = p1.xi() - p2.xi();
    scalar_type const length   = diff.norm();
    Eigen::Vector3d const n    = diff / length;
    auto const C               = length - d_;

    scalar_type const weighted_sum_of_gradients = w1 + w2;
    scalar_type const alpha_tilde               = alpha_ / (dt * dt);
    scalar_type const delta_lagrange =
        -(C + alpha_tilde * lagrange_) / (weighted_sum_of_gradients + alpha_tilde);

    lagrange_ += delta_lagrange;
    p1.xi() += w1 * n * delta_lagrange;
    p2.xi() += w2 * -n * delta_lagrange;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs