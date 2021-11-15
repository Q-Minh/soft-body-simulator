#ifndef SBS_PHYSICS_XPBD_HYBRID_MESH_MESHLESS_MLS_SURFACE_CONTACT_HANDLER_H
#define SBS_PHYSICS_XPBD_HYBRID_MESH_MESHLESS_MLS_SURFACE_CONTACT_HANDLER_H

#include "sbs/physics/collision/contact.h"

namespace sbs {
namespace physics {
namespace xpbd {

// Forward declares
class simulation_t;

class hybrid_mesh_meshless_mls_surface_contact_handler_t : public collision::contact_handler_t
{
  public:
    hybrid_mesh_meshless_mls_surface_contact_handler_t(simulation_t& simulation);

    virtual void on_cd_starting() override;
    virtual void on_cd_ending() override;
    virtual void handle(collision::contact_t const& contact) override;

  private:
    simulation_t& simulation_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_HYBRID_MESH_MESHLESS_MLS_SURFACE_CONTACT_HANDLER_H