#include <sbs/physics/particle.h>

namespace sbs {
namespace physics {

particle_t::particle_t(position_type const& p)
    : x0_(p), xi_(p), xn_(p), x_(p), v_(0., 0., 0.), f_(0., 0., 0.), m_(1.)
{
}

particle_t::position_type const& particle_t::x0() const
{
    return x0_;
}
particle_t::position_type const& particle_t::x() const
{
    return x_;
}
particle_t::position_type const& particle_t::xi() const
{
    return xi_;
}
particle_t::position_type const& particle_t::xn() const
{
    return xn_;
}
particle_t::velocity_type const& particle_t::v() const
{
    return v_;
}
particle_t::force_type const& particle_t::f() const
{
    return f_;
}
scalar_type const& particle_t::mass() const
{
    return m_;
}
scalar_type particle_t::invmass() const
{
    scalar_type constexpr zero{0.};
    scalar_type constexpr one{1.};
    return m_ > zero ? one / m_ : zero;
}
bool particle_t::fixed() const
{
    scalar_type constexpr zero{0.};
    return m_ == zero;
}
particle_t::acceleration_type particle_t::a() const
{
    return f_ * invmass();
}

particle_t::position_type& particle_t::x0()
{
    return x0_;
}
particle_t::position_type& particle_t::x()
{
    return x_;
}
particle_t::position_type& particle_t::xi()
{
    return xi_;
}
particle_t::position_type& particle_t::xn()
{
    return xn_;
}
particle_t::velocity_type& particle_t::v()
{
    return v_;
}
particle_t::force_type& particle_t::f()
{
    return f_;
}
scalar_type& particle_t::mass()
{
    return m_;
}

} // namespace physics
} // namespace sbs