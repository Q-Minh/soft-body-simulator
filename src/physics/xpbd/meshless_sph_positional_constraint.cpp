#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_node.h>
#include <sbs/physics/mechanics/meshless_sph_body.h>
#include <sbs/physics/mechanics/meshless_sph_surface.h>
#include <sbs/physics/simulation.h>
#include <sbs/physics/xpbd/meshless_sph_positional_constraint.h>

namespace sbs {
namespace physics {
namespace xpbd {

meshless_sph_positional_constraint_t::meshless_sph_positional_constraint_t(
    scalar_type alpha,
    scalar_type beta,
    index_type bi,
    mechanics::meshless_sph_body_t const* b,
    mechanics::meshless_sph_surface_vertex_t const& vk,
    Eigen::Vector3d const& vkp,
    scalar_type mu,
    Eigen::Vector3d const& pi)
    : constraint_t(alpha, beta), bi_(bi), b_(b), vk_(vk), vkp_(vkp), mu_(mu), pi_(pi)
{
}

void meshless_sph_positional_constraint_t::project_positions(
    simulation_t& simulation,
    scalar_type dt)
{
    auto& particles = simulation.particles()[bi_];

    auto const& neighbours = vk_.neighbours();
    auto const& Xkjs       = vk_.Xkjs();
    auto const& Wkjs       = vk_.Wkjs();
    auto const& Vjs        = vk_.Vjs();
    auto const sk          = vk_.sk();

    Eigen::Vector3d const diff = vkp_ - pi_;
    scalar_type const C        = scalar_type{0.5} * mu_ * diff.squaredNorm();

    Eigen::Matrix3d Fi{};

    std::vector<Eigen::Vector3d> gradC{};
    gradC.reserve(neighbours.size());
    scalar_type weighted_sum_of_gradients{0.};
    scalar_type gradC_dot_displacement{0.};
    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        index_type const j         = neighbours[a];
        scalar_type const& Vj      = Vjs[a];
        scalar_type const& Wkj     = Wkjs[a];
        Eigen::Vector3d const grad = mu_ * diff * sk * Vj * Wkj;
        particle_t const& p        = particles[j];
        weighted_sum_of_gradients += p.invmass() * grad.squaredNorm();
        gradC.push_back(grad);
        gradC_dot_displacement += gradC[a].dot(p.xi() - p.xn());
    }

    scalar_type const dt2         = dt * dt;
    scalar_type const alpha_tilde = alpha() / dt2;
    scalar_type const beta_tilde  = beta() * dt2;
    scalar_type const gamma       = alpha_tilde * beta_tilde / dt;

    scalar_type const delta_lagrange_num =
        -(C + alpha_tilde * lagrange_) - gamma * gradC_dot_displacement;
    scalar_type const delta_lagrange_den = (1. + gamma) * weighted_sum_of_gradients + alpha_tilde;
    scalar_type const delta_lagrange     = delta_lagrange_num / delta_lagrange_den;

    vkp_.setZero();
    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        index_type const j = neighbours[a];
        particle_t& p      = particles[j];
        p.xi() += p.invmass() * gradC[a] * delta_lagrange;

        scalar_type const& Vj      = Vjs[a];
        Eigen::Vector3d const& Xkj = Xkjs[a];
        scalar_type const& Wkj     = Wkjs[a];
        Eigen::Matrix3d const Fj   = b_->nodes()[j].Fi();

        // Update vertex position
        vkp_ += sk * Vj * (Fj * Xkj + p.xi()) * Wkj;
    }
}

} // namespace xpbd
} // namespace physics
} // namespace sbs
