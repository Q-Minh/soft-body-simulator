#include "sbs/physics/xpbd/constraint.h"

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

void constraint_t::set_as_active()
{
    on_ = true;
}

void constraint_t::set_as_inactive()
{
    on_ = false;
}

bool constraint_t::is_active()
{
    return on_;
}

void constraint_t::toggle()
{
    on_ = !on_;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs
