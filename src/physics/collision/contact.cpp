#include <sbs/physics/collision/contact.h>

namespace sbs {
namespace physics {
namespace collision {

contact_t::type_t contact_t::type() const
{
    return type_;
}
Eigen::Vector3d const& contact_t::point() const
{
    return point_;
}
Eigen::Vector3d const& contact_t::normal() const
{
    return normal_;
}
index_type contact_t::b1() const
{
    return bodies_[0u];
}
index_type contact_t::b2() const
{
    return bodies_[1u];
}

index_type surface_mesh_particle_to_sdf_contact_t::vi() const
{
    return vi_;
}

index_type& surface_mesh_particle_to_sdf_contact_t::vi()
{
    return vi_;
}

} // namespace collision
} // namespace physics
} // namespace sbs