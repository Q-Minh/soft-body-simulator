#ifndef SBS_PHYSICS_COLLISION_CD_SYSTEM_H
#define SBS_PHYSICS_COLLISION_CD_SYSTEM_H

#include <memory>
#include <sbs/physics/collision/contact.h>
#include <utility>
#include <vector>

namespace sbs {
namespace physics {
namespace xpbd {
class simulation_t;
} // namespace xpbd

namespace collision {

// Forward declares
class collision_model_t;

class cd_system_t
{
  public:
    using intersection_pair = std::pair<collision_model_t*, collision_model_t*>;

    cd_system_t(std::vector<collision_model_t*> const& collision_objects);

    cd_system_t(cd_system_t const& other) = default;
    cd_system_t(cd_system_t&& other)      = default;

    cd_system_t& operator=(cd_system_t const& other) = default;
    cd_system_t& operator=(cd_system_t&& other) = default;

    virtual void execute()                                    = 0;
    virtual void update(xpbd::simulation_t const& simulation) = 0;

    std::vector<collision_model_t*> const& collision_objects() const;

    std::unique_ptr<contact_handler_t> const& contact_handler() const;
    std::unique_ptr<contact_handler_t>& contact_handler();
    void use_contact_handler(std::unique_ptr<contact_handler_t> contact_handler);

  protected:
    std::vector<collision_model_t*>& collision_objects();

  private:
    std::vector<collision_model_t*> collision_objects_;
    std::unique_ptr<contact_handler_t> contact_handler_{nullptr};
};

} // namespace collision
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_CD_SYSTEM_H