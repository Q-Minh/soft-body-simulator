#include "physics/xpbd/collision_constraint.h"

#include <Eigen/Geometry>

namespace sbs {
namespace physics {
namespace xpbd {

collision_constraint_t::collision_constraint_t(
    //scalar_type const alpha,
    index_pair_type penetrating_vertex,
    position_type const& q,
    normal_type const& n)
    : constraint_t(0.),
      p_(penetrating_vertex.first),
      b1_(penetrating_vertex.second),
      q_(q),
      n_(n)
{
}

void collision_constraint_t::project(
    std::vector<positions_type>& positions,
    std::vector<masses_type> const& masses,
    scalar_type& lagrange_multiplier,
    scalar_type const dt) const
{
    Eigen::Vector3d const p = positions[b1_].col(p_);
    scalar_type const m     = masses[b1_](p_);
    scalar_type const w     = 1. / m;

    scalar_type const C = evaluate(p);

    /**
     * Inequality constraint C >= 0 is already satisfied
     */
    if (C >= 0.)
        return;

    /**
     * grad(C) = n
     * ||n||^2 = 1,
     * so w*||n||^2 = w
     */
    scalar_type const alpha_tilde    = alpha_ / (dt * dt);
    scalar_type const delta_lagrange = -(C + alpha_tilde * lagrange_multiplier) / (w + alpha_tilde);

    lagrange_multiplier += delta_lagrange;
    positions[b1_].col(p_) += w * n_ * delta_lagrange;
}

collision_constraint_t::scalar_type collision_constraint_t::evaluate(position_type const& p) const
{
    Eigen::Vector3d const qp = p - q_;
    return qp.dot(n_);
}

} // namespace xpbd
} // namespace physics
} // namespace sbs