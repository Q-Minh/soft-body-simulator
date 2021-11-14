#include "sbs/physics/body/body.h"

namespace sbs {
namespace physics {
namespace body {

index_type body_t::id() const
{
    return id_;
}

index_type& body_t::id()
{
    return id_;
}

} // namespace body
} // namespace physics
} // namespace sbs
