#include "sbs/physics/xpbd/contact_handler.h"

#include "sbs/physics/body/body.h"
#include "sbs/physics/tetrahedral_body.h"
#include "sbs/physics/tetrahedral_mesh_boundary.h"
#include "sbs/physics/xpbd/collision_constraint.h"
#include "sbs/physics/xpbd/simulation.h"

namespace sbs {
namespace physics {
namespace xpbd {

contact_handler_t::contact_handler_t(simulation_t& simulation) : simulation_(simulation) {}

void contact_handler_t::on_cd_starting()
{
    // no-op
}

void contact_handler_t::on_cd_ending()
{
    // no-op
}

void contact_handler_t::handle(collision::contact_t const& contact)
{
    auto const contact_type = contact.type();
    if (contact_type == collision::contact_t::type_t::surface_particle_to_sdf)
    {
        auto const& surface_vertex_to_sdf_contact =
            reinterpret_cast<collision::surface_mesh_particle_to_sdf_contact_t const&>(contact);

        if (on_mesh_vertex_to_sdf_contact)
        {
            on_mesh_vertex_to_sdf_contact(surface_vertex_to_sdf_contact);
        }
    }
}

} // namespace xpbd
} // namespace physics
} // namespace sbs