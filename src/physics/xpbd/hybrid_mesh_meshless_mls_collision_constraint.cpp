#include <optional>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_mls_body.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_mls_node.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_mls_surface.h>
#include <sbs/physics/simulation.h>
#include <sbs/physics/xpbd/hybrid_mesh_meshless_mls_collision_constraint.h>

namespace sbs {
namespace physics {
namespace xpbd {

hybrid_mesh_meshless_mls_collision_constraint_t::hybrid_mesh_meshless_mls_collision_constraint_t(
    scalar_type alpha,
    scalar_type beta,
    index_type const bi,
    mechanics::hybrid_mesh_meshless_mls_body_t const* b,
    mechanics::hybrid_mesh_meshless_mls_surface_vertex_t const& vk,
    Eigen::Vector3d const& vi,
    Eigen::Vector3d const& c,
    Eigen::Vector3d const& n)
    : constraint_t(alpha, beta), bi_(bi), b_(b), vk_(vk), vi_(vi), c_(c), n_(n)
{
}

void hybrid_mesh_meshless_mls_collision_constraint_t::project_positions(
    simulation_t& simulation,
    scalar_type dt)
{
    auto& particles = simulation.particles()[bi_];

    auto const mesh_particles_index_offset     = b_->get_mesh_particles_index_offset();
    auto const meshless_particles_index_offset = b_->get_meshless_particles_index_offset();

    auto const& neighbours = vk_.neighbours();

    scalar_type const C = (vi_ - c_).dot(n_);

    if (C >= 0.)
    {
        return;
    }

    scalar_type weighted_sum_of_gradients{0.};
    scalar_type gradC_dot_displacement{0.};

    std::vector<Eigen::Vector3d> gradC_meshless{};
    gradC_meshless.reserve(neighbours.size());
    std::vector<scalar_type> const& phi_js = vk_.phi_js();
    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        index_type const j         = neighbours[a];
        std::size_t const pidx     = static_cast<std::size_t>(meshless_particles_index_offset) + j;
        particle_t const& pj       = particles[pidx];
        scalar_type const phi_j    = phi_js[a];
        Eigen::Vector3d const grad = phi_j * n_;
        weighted_sum_of_gradients += pj.invmass() * grad.squaredNorm();
        gradC_dot_displacement += grad.dot(pj.xi() - pj.xn());
        gradC_meshless.push_back(grad);
    }

    // Compute mesh shape function contributions if the surface vertex
    // is influenced by mesh shape functions
    std::array<std::optional<Eigen::Vector3d>, 4u> gradC_mesh{};
    if (vk_.is_in_tetrahedron())
    {
        std::array<std::optional<scalar_type>, 4u> const& mesh_phi_js = vk_.mesh_phi_js();
        topology::tetrahedron_t const& t = b_->topology().tetrahedron(vk_.ti());
        auto const& vis        = t.vertex_indices();
        for (std::uint8_t v = 0u; v < 4u; ++v)
        {
            if (!mesh_phi_js[v].has_value())
                continue;

            index_type const vi           = vis[v];
            scalar_type const& mesh_phi_j = mesh_phi_js[v].value();
            Eigen::Vector3d const grad    = mesh_phi_j * n_;
            gradC_mesh[v]                 = grad;
            std::size_t const pidx = static_cast<std::size_t>(mesh_particles_index_offset) + vi;
            particle_t const& pj   = particles[pidx];
            weighted_sum_of_gradients += pj.invmass() * grad.squaredNorm();
            gradC_dot_displacement += grad.dot(pj.xi() - pj.xn());
        }
    }

    scalar_type const dt2         = dt * dt;
    scalar_type const alpha_tilde = alpha() / dt2;
    scalar_type const beta_tilde  = beta() * dt2;
    scalar_type const gamma       = alpha_tilde * beta_tilde / dt;

    scalar_type const delta_lagrange_num =
        -(C + alpha_tilde * lagrange_) - gamma * gradC_dot_displacement;
    scalar_type const delta_lagrange_den = (1. + gamma) * weighted_sum_of_gradients + alpha_tilde;
    scalar_type const delta_lagrange     = delta_lagrange_num / delta_lagrange_den;

    // Perform position projection and also recompute the new surface vertex position
    vi_.setZero();
    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        index_type const j     = neighbours[a];
        std::size_t const pidx = static_cast<std::size_t>(meshless_particles_index_offset) + j;
        particle_t& pj         = particles[pidx];
        pj.xi() += pj.invmass() * gradC_meshless[a] * delta_lagrange;

        scalar_type const& phi_j = phi_js[a];
        vi_ += pj.xi() * phi_j;
    }
    if (vk_.is_in_tetrahedron())
    {
        std::array<std::optional<scalar_type>, 4u> const& mesh_phi_js = vk_.mesh_phi_js();
        topology::tetrahedron_t const& t = b_->topology().tetrahedron(vk_.ti());
        auto const& vis        = t.vertex_indices();
        for (std::uint8_t v = 0u; v < 4u; ++v)
        {
            if (!gradC_mesh[v].has_value())
                continue;

            index_type const vi    = vis[v];
            std::size_t const pidx = static_cast<std::size_t>(mesh_particles_index_offset) + vi;
            particle_t& pj         = particles[pidx];
            Eigen::Vector3d const& grad = gradC_mesh[v].value();
            pj.xi() += pj.invmass() * grad * delta_lagrange;

            scalar_type const& mesh_phi_j = mesh_phi_js[v].value();
            vi_ += pj.xi() * mesh_phi_j;
        }
    }
}

} // namespace xpbd
} // namespace physics
} // namespace sbs
