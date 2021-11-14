#ifndef SBS_PHYSICS_XPBD_CONTACT_HANDLER_H
#define SBS_PHYSICS_XPBD_CONTACT_HANDLER_H

#include "sbs/physics/collision/contact.h"

#include <functional>

namespace sbs {
namespace physics {
namespace xpbd {

class simulation_t;

class contact_handler_t : public collision::contact_handler_t
{
  public:
    contact_handler_t(simulation_t& simulation);

    virtual void on_cd_starting() override;
    virtual void on_cd_ending() override;
    virtual void handle(collision::contact_t const& contact) override;

    std::function<void(collision::surface_mesh_particle_to_sdf_contact_t const& contact)>
        on_mesh_vertex_to_sdf_contact;

  private:
    simulation_t& simulation_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_CONTACT_HANDLER_H
