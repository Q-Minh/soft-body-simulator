#include <sbs/physics/body.h>

namespace sbs {
namespace physics {

index_type body_t::id() const
{
    return id_;
}

index_type& body_t::id()
{
    return id_;
}

} // namespace physics
} // namespace sbs
