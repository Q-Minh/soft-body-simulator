#include "physics/xpbd/constraint.h"

namespace sbs {
namespace physics {
namespace xpbd {

constraint_t::scalar_type const& constraint_t::alpha() const
{
    return alpha_;
}

constraint_t::scalar_type& constraint_t::alpha()
{
    return alpha_;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs
