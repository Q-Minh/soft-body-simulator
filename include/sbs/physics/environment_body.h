#ifndef SBS_PHYSICS_ENVIRONMENT_BODY_H
#define SBS_PHYSICS_ENVIRONMENT_BODY_H

#include <sbs/physics/body.h>

namespace sbs {
namespace physics {

class environment_body_t : public body_t
{
  public:
    virtual visual_model_type const& visual_model() const override;
    virtual collision_model_type const& collision_model() const override;
    virtual void update_visual_model(simulation_t const& simulation) override;
    virtual void update_collision_model(simulation_t const& simulation) override;
    virtual void update_physical_model(simulation_t const& simulation) override;

  private:
    
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_ENVIRONMENT_BODY_H