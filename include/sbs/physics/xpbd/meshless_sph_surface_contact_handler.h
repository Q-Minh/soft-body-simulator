#ifndef SBS_PHYSICS_XPBD_MESHLESS_SPH_SURFACE_CONTACT_HANDLER_H
#define SBS_PHYSICS_XPBD_MESHLESS_SPH_SURFACE_CONTACT_HANDLER_H

#include <sbs/physics/collision/contact.h>

namespace sbs {
namespace physics {

class simulation_t;

namespace xpbd {

class meshless_sph_surface_contact_handler_t : public collision::contact_handler_t
{
  public:
    meshless_sph_surface_contact_handler_t(simulation_t& simulation);

    virtual void on_cd_starting() override;
    virtual void on_cd_ending() override;
    virtual void handle(collision::contact_t const& contact) override;

  private:
    simulation_t& simulation_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_MESHLESS_SPH_SURFACE_CONTACT_HANDLER_H
