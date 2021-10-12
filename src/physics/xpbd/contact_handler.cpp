#include <sbs/physics/body.h>
#include <sbs/physics/simulation.h>
#include <sbs/physics/tetrahedral_body.h>
#include <sbs/physics/tetrahedral_mesh_boundary.h>
#include <sbs/physics/xpbd/collision_constraint.h>
#include <sbs/physics/xpbd/contact_handler.h>

namespace sbs {
namespace physics {
namespace xpbd {

void contact_handler_t::handle(collision::contact_t const& contact)
{
    auto const contact_type = contact.type();
    if (contact_type == collision::contact_t::type_t::surface_particle_to_sdf)
    {
        collision::surface_mesh_particle_contact_t const& surface_mesh_contact =
            reinterpret_cast<collision::surface_mesh_particle_contact_t const&>(contact);

        body_t const& b1 = *simulation_.bodies().at(contact.b1());
        body_t const& b2 = *simulation_.bodies().at(contact.b2());

        tetrahedral_body_t const* tet_mesh = dynamic_cast<tetrahedral_body_t const*>(&b1);
        if (!tet_mesh)
        {
            return;
        }

        auto const& visual_model = b1.visual_model();
        tetrahedral_mesh_boundary_t const* mesh_boundary =
            dynamic_cast<tetrahedral_mesh_boundary_t const*>(&visual_model);
        if (!mesh_boundary)
        {
            return;
        }

        index_type const particle_index =
            mesh_boundary->from_surface_vertex(surface_mesh_contact.vi());

        auto collision_constraint = std::make_unique<xpbd::collision_constraint_t>(
            simulation_.simulation_parameters().compliance,
            simulation_.simulation_parameters().collision_damping,
            simulation_,
            b1.id(),
            particle_index,
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