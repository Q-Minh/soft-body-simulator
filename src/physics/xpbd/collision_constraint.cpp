#include "sbs/physics/xpbd/collision_constraint.h"

#include "sbs/common/primitive.h"
#include "sbs/physics/xpbd/mesh.h"

#include <Eigen/Geometry>

namespace sbs {
namespace physics {
namespace xpbd {

collision_constraint_t::collision_constraint_t(
    scalar_type alpha,
    body_ptr_type b,
    index_type vi,
    common::triangle_t const& t,
    normal_type const& n)
    : constraint_t{alpha}, b_(b), vi_(vi), qs_(), n_(n)
{
    auto const& p = b_->vertices().at(vi_).position();
    qs_           = common::closest_point(p, t);
}

void collision_constraint_t::project(
    std::vector<std::shared_ptr<tetrahedral_mesh_t>> const& bodies,
    scalar_type& lagrange_multiplier,
    scalar_type const dt) const
{
    physics::vertex_t& v    = b_->vertices().at(vi_);
    Eigen::Vector3d const p = v.position();
    scalar_type const w     = 1. / v.mass();

    scalar_type const C = evaluate(p);

    /**
     * Inequality constraint C >= 0 is already satisfied
     */
    if (C >= 0.)
        return;

    scalar_type const alpha_tilde    = alpha() / (dt * dt);
    scalar_type const delta_lagrange = -(C + alpha_tilde * lagrange_multiplier) / (w + alpha_tilde);

    lagrange_multiplier += delta_lagrange;

    Eigen::Vector3d const correction_direction = n_; // try some other direction?

    /**
     * grad(C) = n
     * ||n||^2 = 1,
     * so w*||n||^2 = w
     */
    v.position() += w * correction_direction * delta_lagrange;
}

collision_constraint_t::scalar_type collision_constraint_t::evaluate(position_type const& p) const
{
    Eigen::Vector3d const qp = p - qs_;
    double const C           = qp.dot(n_);
    return C;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs
