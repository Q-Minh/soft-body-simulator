#include <sbs/physics/collision/cd_system.h>

namespace sbs {
namespace physics {
namespace collision {

cd_system_t::cd_system_t(std::vector<collision_model_t*> const& collision_objects)
    : collision_objects_(collision_objects)
{
}
std::vector<collision_model_t*> const& cd_system_t::collision_objects() const
{
    return collision_objects_;
}

std::unique_ptr<contact_handler_t> const& cd_system_t::contact_handler() const
{
    return contact_handler_;
}
std::unique_ptr<contact_handler_t>& cd_system_t::contact_handler()
{
    return contact_handler_;
}
void cd_system_t::use_contact_handler(std::unique_ptr<contact_handler_t> contact_handler)
{
    contact_handler_ = std::move(contact_handler);
}
std::vector<collision_model_t*>& cd_system_t::collision_objects()
{
    return collision_objects_;
}

} // namespace collision
} // namespace physics
} // namespace sbs