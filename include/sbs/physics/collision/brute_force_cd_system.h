#ifndef SBS_PHYSICS_COLLISION_BRUTE_FORCE_CD_SYSTEM_H
#define SBS_PHYSICS_COLLISION_BRUTE_FORCE_CD_SYSTEM_H

#include <sbs/physics/collision/cd_system.h>

namespace sbs {
namespace physics {
namespace collision {

class brute_force_cd_system_t : public cd_system_t
{
  public:
    brute_force_cd_system_t(std::vector<collision_model_t*> const& collision_objects);

    virtual void execute() override;
    virtual void update(xpbd::simulation_t const& simulation) override;
};

} // namespace collision
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_BRUTE_FORCE_CD_SYSTEM_H