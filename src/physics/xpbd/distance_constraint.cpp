#include "sbs/physics/xpbd/distance_constraint.h"

#include "sbs/physics/xpbd/mesh.h"

namespace sbs {
namespace physics {
namespace xpbd {

distance_constraint_t::distance_constraint_t(
    scalar_type const alpha,
    position_key_type const& vb1,
    position_key_type const& vb2)
    : constraint_t(alpha), b1_(vb1.first), v1_(vb1.second), b2_(vb2.first), v2_(vb2.second), d_()
{
    Eigen::Vector3d const p1 = b1_->vertices().at(v1_).position();
    Eigen::Vector3d const p2 = b2_->vertices().at(v2_).position();
    d_                       = (p2 - p1).norm();
}

void distance_constraint_t::project(
    std::vector<std::shared_ptr<xpbd::tetrahedral_mesh_t>> const& bodies,
    scalar_type& lagrange_multiplier,
    scalar_type const dt) const
{
    physics::vertex_t& v1 = b1_->vertices().at(v1_);
    physics::vertex_t& v2 = b2_->vertices().at(v2_);

    Eigen::Vector3d const p1 = v1.position();
    Eigen::Vector3d const p2 = v2.position();

    double const m1 = v1.mass();
    double const m2 = v2.mass();

    double const w1 = v1.fixed() ? 0.0 : 1. / m1;
    double const w2 = v2.fixed() ? 0.0 : 1. / m2;

    auto const n = (p1 - p2).normalized();
    auto const C = evaluate(p1, p2);

    Eigen::Vector3d const diff = p2 - p1;
    scalar_type const length   = diff.norm();

    // <n,n> = 1 and <-n,-n> = 1
    scalar_type const weighted_sum_of_gradients = w1 + w2;
    scalar_type const alpha_tilde               = alpha() / (dt * dt);
    scalar_type const delta_lagrange =
        -(C + alpha_tilde * lagrange_multiplier) / (weighted_sum_of_gradients + alpha_tilde);

    lagrange_multiplier += delta_lagrange;
    v1.position() += w1 * n * delta_lagrange;
    v2.position() += w2 * -n * delta_lagrange;
}

distance_constraint_t::scalar_type
distance_constraint_t::evaluate(Eigen::Vector3d const& p1, Eigen::Vector3d const& p2) const
{
    return (p2 - p1).norm() - d_;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs