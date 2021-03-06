#include <sbs/physics/collision/collision_model.h>

namespace sbs {
namespace physics {
namespace collision {

collision_model_t::volume_type const& collision_model_t::volume() const
{
    return englobing_volume_;
}

collision_model_t::volume_type& collision_model_t::volume()
{
    return englobing_volume_;
}

index_type collision_model_t::id() const
{
    return id_;
}

index_type& collision_model_t::id()
{
    return id_;
}

} // namespace collision
} // namespace physics
} // namespace sbs
