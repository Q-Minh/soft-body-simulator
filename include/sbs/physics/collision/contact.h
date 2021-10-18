#ifndef SBS_PHYSICS_COLLISION_CONTACT_H
#define SBS_PHYSICS_COLLISION_CONTACT_H

#include <Eigen/Core>
#include <sbs/aliases.h>

namespace sbs {
namespace physics {
namespace collision {

class contact_t
{
  public:
    enum class type_t { surface_particle_to_sdf };

    contact_t(
        type_t contact_type,
        index_type const body1,
        index_type const body2,
        Eigen::Vector3d const& contact_point,
        Eigen::Vector3d const& contact_normal)
        : type_(contact_type), bodies_{body1, body2}, point_(contact_point), normal_(contact_normal)
    {
    }

    type_t type() const;
    Eigen::Vector3d const& point() const;
    Eigen::Vector3d const& normal() const;
    index_type b1() const;
    index_type b2() const;

  private:
    type_t type_;
    index_type bodies_[2];
    Eigen::Vector3d point_;
    Eigen::Vector3d normal_;
};

class surface_mesh_particle_to_sdf_contact_t : public contact_t
{
  public:
    surface_mesh_particle_to_sdf_contact_t(
        contact_t::type_t contact_type,
        index_type const body1,
        index_type const body2,
        Eigen::Vector3d const& contact_point,
        Eigen::Vector3d const& contact_normal,
        index_type const vi)
        : contact_t(contact_type, body1, body2, contact_point, contact_normal), vi_(vi)
    {
    }

    index_type vi() const;
    index_type& vi();

  private:
    index_type vi_;
};

class contact_handler_t
{
  public:
    virtual void on_cd_starting()                 = 0;
    virtual void on_cd_ending()                   = 0;
    virtual void handle(contact_t const& contact) = 0;
};

} // namespace collision
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_CONTACT_H