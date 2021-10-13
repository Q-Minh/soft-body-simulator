#include <sbs/physics/collision/brute_force_cd_system.h>
#include <sbs/physics/collision/collision_model.h>

namespace sbs {
namespace physics {
namespace collision {

brute_force_cd_system_t::brute_force_cd_system_t(
    std::vector<collision_model_t*> const& collision_objects)
    : cd_system_t(collision_objects)
{
}

void brute_force_cd_system_t::execute()
{
    std::vector<collision_model_t*>& objects = collision_objects();
    for (std::size_t i = 0u; i < objects.size(); ++i)
    {
        collision_model_t* b1 = objects[i];
        for (std::size_t j = 0u; j < objects.size(); ++j)
        {
            collision_model_t* b2 = objects[j];
            b1->collide(*b2, *contact_handler());
        }
    }
}

void brute_force_cd_system_t::update(simulation_t const& simulation)
{
    // no-op
}

} // namespace collision
} // namespace physics
} // namespace sbs