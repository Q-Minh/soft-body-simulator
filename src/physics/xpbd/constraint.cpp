#include "sbs/physics/xpbd/constraint.h"

#include "sbs/physics/xpbd/simulation.h"

namespace sbs {
namespace physics {
namespace xpbd {

constraint_t::constraint_t(scalar_type alpha, scalar_type beta)
    : alpha_(alpha), beta_(beta), lagrange_(0.), js_(), bis_()
{
}

constraint_t::constraint_t(
    scalar_type alpha,
    scalar_type beta,
    std::vector<index_type> const& js,
    std::vector<index_type> const& bis)
    : alpha_(alpha), beta_(beta), lagrange_(0.), js_(js), bis_(bis)
{
    assert(js_.size() == bis_.size());
}

void constraint_t::prepare_for_projection(simulation_t& simulation)
{
    lagrange_ = 0.;
    prepare_for_projection_impl(simulation);
}

void constraint_t::project_positions_with_dampling(
    simulation_t& simulation,
    scalar_type const C,
    std::vector<Eigen::Vector3d> const& gradC,
    scalar_type dt)
{
    scalar_type weighted_sum_of_gradients{0.};
    scalar_type gradC_dot_displacement{0.};

    assert(js_.size() == gradC.size());
    assert(bis_.size() == js_.size());

    auto& particles = simulation.particles();
    for (auto i = 0u; i < js_.size(); ++i)
    {
        index_type const j   = js_[i];
        index_type const b   = bis_[i];
        particle_t const& pj = particles[b][j];
        weighted_sum_of_gradients += pj.invmass() * gradC[i].squaredNorm();
        gradC_dot_displacement += gradC[i].dot(pj.xi() - pj.xn());
    }

    scalar_type const dt2         = dt * dt;
    scalar_type const alpha_tilde = alpha() / dt2;
    scalar_type const beta_tilde  = beta() * dt2;
    scalar_type const gamma       = alpha_tilde * beta_tilde / dt;

    scalar_type const delta_lagrange_num =
        -(C + alpha_tilde * lagrange_) - gamma * gradC_dot_displacement;
    scalar_type const delta_lagrange_den = (1. + gamma) * (weighted_sum_of_gradients) + alpha_tilde;
    scalar_type const delta_lagrange     = delta_lagrange_num / delta_lagrange_den;

    lagrange_ += delta_lagrange;

    for (auto i = 0u; i < js_.size(); ++i)
    {
        index_type const j = js_[i];
        index_type const b = bis_[i];
        particle_t& pj     = particles[b][j];
        pj.xi() += pj.invmass() * gradC[i] * delta_lagrange;
    }
}

scalar_type constraint_t::alpha() const
{
    return alpha_;
}
scalar_type constraint_t::beta() const
{
    return beta_;
}
scalar_type constraint_t::lambda() const
{
    return lagrange_;
}

scalar_type constraint_t::compliance() const
{
    return alpha_;
}
scalar_type constraint_t::damping() const
{
    return beta_;
}

std::vector<index_type> const& constraint_t::js() const
{
    return js_;
}

std::vector<index_type> const& constraint_t::bs() const
{
    return bis_;
}
std::vector<index_type>& constraint_t::js()
{
    return js_;
}

std::vector<index_type>& constraint_t::bs()
{
    return bis_;
}

void constraint_t::set_indices(std::vector<index_type> const& js)
{
    js_ = js;
}

void constraint_t::set_bodies(std::vector<index_type> const& bis)
{
    bis_ = bis;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs