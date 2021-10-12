#include <sbs/physics/constraint.h>
#include <sbs/physics/simulation.h>

namespace sbs {
namespace physics {

constraint_t::constraint_t(scalar_type alpha, scalar_type beta)
    : alpha_(alpha), beta_(beta), lagrange_(0.)
{
}

void constraint_t::prepare_for_projection(simulation_t& simulation)
{
    lagrange_ = 0.;
    prepare_for_projection_impl(simulation);
}

scalar_type constraint_t::alpha() const
{
    return alpha_;
}
scalar_type constraint_t::beta() const
{
    return beta_;
}
scalar_type constraint_t::lambda() const
{
    return lagrange_;
}

scalar_type constraint_t::compliance() const
{
    return alpha_;
}
scalar_type constraint_t::damping() const
{
    return beta_;
}

} // namespace physics
} // namespace sbs