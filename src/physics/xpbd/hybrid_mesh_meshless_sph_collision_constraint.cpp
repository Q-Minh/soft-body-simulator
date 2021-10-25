#include <optional>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_body.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_node.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_surface.h>
#include <sbs/physics/simulation.h>
#include <sbs/physics/xpbd/hybrid_mesh_meshless_sph_collision_constraint.h>

namespace sbs {
namespace physics {
namespace xpbd {

sbs::physics::xpbd::hybrid_mesh_meshless_mls_collision_constraint_t::
    hybrid_mesh_meshless_mls_collision_constraint_t(
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

    auto const mesh_particles_index_offset     = 0u;
    auto const meshless_particles_index_offset = b_->mesh_node_count();

    auto const& neighbours = vk_.neighbours();
    auto const& Xkjs       = vk_.Xkjs();
    auto const& Wkjs       = vk_.Wkjs();
    auto const& Vjs        = vk_.Vjs();
    auto const sk          = vk_.sk();

    scalar_type const C = (vi_ - c_).dot(n_);

    if (C >= 0.)
    {
        return;
    }

    std::vector<Eigen::Vector3d> gradC_meshless{};
    gradC_meshless.reserve(neighbours.size());
    scalar_type weighted_sum_of_gradients{0.};
    scalar_type gradC_dot_displacement{0.};
    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        index_type const j         = neighbours[a];
        scalar_type const& Vj      = Vjs[a];
        scalar_type const& Wkj     = Wkjs[a];
        Eigen::Vector3d const grad = sk * Vj * Wkj * n_;
        particle_t const& p        = particles[meshless_particles_index_offset + j];
        weighted_sum_of_gradients += p.invmass() * grad.squaredNorm();
        gradC_meshless.push_back(grad);
        gradC_dot_displacement += gradC_meshless[a].dot(p.xi() - p.xn());
    }

    std::array<std::optional<Eigen::Vector3d>, 4u> gradC_mesh{};
    index_type const ti = vk_.ti();
    bool const is_vertex_in_boundary_mesh =
        ti != std::numeric_limits<index_type>::max() && b_->is_boundary_mesh_tetrahedron(ti);
    if (is_vertex_in_boundary_mesh)
    {
        // Compute contribution of interpolation over mesh nodes
        tetrahedron_t const& t = b_->topology().tetrahedron(ti);
        for (std::uint8_t v = 0u; v < 4u; ++v)
        {
            index_type const vi = t.vertex_indices()[v];
            if (b_->is_boundary_mesh_vertex(vi))
                continue;

            particle_t const& pi = particles[mesh_particles_index_offset + vi];

            auto const& xi          = b_->x()[vi];
            auto const phi          = b_->phi_i(ti, v);
            auto const& Xi          = vk_.x0();
            scalar_type const phi_i = phi.dot(Eigen::Vector4d{1., Xi.x(), Xi.y(), Xi.z()});

            gradC_mesh[v] = phi_i * n_;
            gradC_dot_displacement += gradC_mesh[v].value().dot(pi.xi() - pi.xn());
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

    vi_.setZero();
    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        index_type const j = neighbours[a];
        particle_t& p      = particles[meshless_particles_index_offset + j];
        p.xi() += p.invmass() * gradC_meshless[a] * delta_lagrange;

        scalar_type const& Vj      = Vjs[a];
        Eigen::Vector3d const& Xkj = Xkjs[a];
        scalar_type const& Wkj     = Wkjs[a];
        Eigen::Matrix3d const Fj   = b_->meshless_nodes()[j].Fi();

        // Update vertex position
        vi_ += sk * Vj * (Fj * Xkj + p.xi()) * Wkj;
    }
    if (is_vertex_in_boundary_mesh)
    {
        tetrahedron_t const& t = b_->topology().tetrahedron(ti);
        for (std::uint8_t v = 0u; v < 4u; ++v)
        {
            // if gradC_mesh[v] has a value, it means this mesh node has a shape function
            // contributing to the surface vertex position
            if (gradC_mesh[v].has_value())
            {
                index_type const vi = t.vertex_indices()[v];
                particle_t& p       = particles[mesh_particles_index_offset + vi];
                p.xi() += p.invmass() * gradC_mesh[v].value() * delta_lagrange;

                auto const& xi          = b_->x()[vi];
                auto const phi          = b_->phi_i(ti, v);
                auto const& Xi          = vk_.x0();
                scalar_type const phi_i = phi.dot(Eigen::Vector4d{1., Xi.x(), Xi.y(), Xi.z()});
                vi_ += p.xi() * phi_i;
            }
        }
    }
}

} // namespace xpbd
} // namespace physics
} // namespace sbs
