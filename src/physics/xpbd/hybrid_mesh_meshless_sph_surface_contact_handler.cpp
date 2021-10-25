#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_body.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_surface.h>
#include <sbs/physics/simulation.h>
#include <sbs/physics/xpbd/hybrid_mesh_meshless_sph_collision_constraint.h>
#include <sbs/physics/xpbd/hybrid_mesh_meshless_sph_surface_contact_handler.h>

namespace sbs {
namespace physics {
namespace xpbd {

sbs::physics::xpbd::hybrid_mesh_meshless_mls_surface_contact_handler_t::
    hybrid_mesh_meshless_mls_surface_contact_handler_t(simulation_t& simulation)
    : simulation_(simulation)
{
}

void sbs::physics::xpbd::hybrid_mesh_meshless_mls_surface_contact_handler_t::on_cd_starting() {}

void sbs::physics::xpbd::hybrid_mesh_meshless_mls_surface_contact_handler_t::on_cd_ending() {}

void sbs::physics::xpbd::hybrid_mesh_meshless_mls_surface_contact_handler_t::handle(
    collision::contact_t const& contact)
{
    auto const contact_type = contact.type();
    if (contact_type == collision::contact_t::type_t::surface_particle_to_sdf)
    {
        collision::surface_mesh_particle_to_sdf_contact_t const& surface_mesh_contact =
            reinterpret_cast<collision::surface_mesh_particle_to_sdf_contact_t const&>(contact);

        body_t const& b1 = *simulation_.bodies().at(contact.b1());
        body_t const& b2 = *simulation_.bodies().at(contact.b2());

        auto const& visual_model = b1.visual_model();
        mechanics::hybrid_mesh_meshless_sph_surface_t const* surface_mesh =
            dynamic_cast<mechanics::hybrid_mesh_meshless_sph_surface_t const*>(&visual_model);
        if (!surface_mesh)
        {
            return;
        }

        auto const vi = surface_mesh_contact.vi();
        mechanics::hybrid_mesh_meshless_mls_surface_vertex_t const& vk =
            surface_mesh->embedded_vertices()[vi];
        Eigen::Vector3d const vip = surface_mesh->vertex(surface_mesh_contact.vi()).position;

        auto collision_constraint =
            std::make_unique<xpbd::hybrid_mesh_meshless_mls_collision_constraint_t>(
                simulation_.simulation_parameters().collision_compliance,
                simulation_.simulation_parameters().collision_damping,
                b1.id(),
                surface_mesh->mechanical_model(),
                vk,
                vip,
                contact.point(),
                contact.normal());
        std::vector<std::unique_ptr<physics::constraint_t>>& collision_constraints =
            simulation_.collision_constraints();
        collision_constraints.push_back(std::move(collision_constraint));
    }
}

} // namespace xpbd
} // namespace physics
} // namespace sbs
