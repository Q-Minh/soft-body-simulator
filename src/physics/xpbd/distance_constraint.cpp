#include "physics/xpbd/distance_constraint.h"

#include "common/node.h"

#include <Eigen/Core>

namespace sbs {
namespace physics {
namespace xpbd {

distance_constraint_t::distance_constraint_t(
    scalar_type const alpha,
    positions_type const& positions,
    index_pair_type const& vb1,
    index_pair_type const& vb2)
    : constraint_t(alpha), v1_(vb1.first), v2_(vb2.first), b1_(vb1.second), b2_(vb2.second), d_()
{
    Eigen::Vector3d const p1 = positions.col(v1_);
    Eigen::Vector3d const p2 = positions.col(v2_);
    d_                       = (p2 - p1).norm();
}

void distance_constraint_t::project(
    std::vector<positions_type>& positions,
    std::vector<masses_type> const& masses,
    scalar_type& lagrange_multiplier,
    scalar_type const dt) const
{
    Eigen::Vector3d const p1 = positions[b1_].col(v1_);
    Eigen::Vector3d const p2 = positions[b2_].col(v2_);

    double const m1 = masses[b1_](v1_);
    double const m2 = masses[b2_](v2_);

    double const w1 = 1. / m1;
    double const w2 = 1. / m2;

    auto const n = (p1 - p2).normalized();
    auto const C = evaluate(p1, p2);

    Eigen::Vector3d const diff = p2 - p1;
    scalar_type const length   = diff.norm();

    // <n,n> = 1 and <-n,-n> = 1
    scalar_type const weighted_sum_of_gradients = w1 + w2;
    scalar_type const alpha_tilde               = alpha_ / (dt * dt);
    scalar_type const delta_lagrange =
        -(C + alpha_tilde * lagrange_multiplier) / (weighted_sum_of_gradients + alpha_tilde);

    lagrange_multiplier += delta_lagrange;
    positions[b1_].col(v1_) += w1 * n * delta_lagrange;
    positions[b2_].col(v2_) += w2 * -n * delta_lagrange;
}

distance_constraint_t::scalar_type
distance_constraint_t::evaluate(Eigen::Vector3d const& p1, Eigen::Vector3d const& p2) const
{
    return (p2 - p1).norm() - d_;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs