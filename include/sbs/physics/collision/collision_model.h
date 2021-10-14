#ifndef SBS_PHYSICS_COLLISION_COLLISION_MODEL_H
#define SBS_PHYSICS_COLLISION_COLLISION_MODEL_H

#include <Eigen/Geometry>
#include <sbs/aliases.h>

namespace sbs {
namespace physics {

// Forward declares
class simulation_t;

namespace collision {

// Forward declares
class contact_handler_t;

enum class model_type_t { bvh, sdf };
enum class bounding_volume_t { aabb, sphere };
enum class primitive_type_t { point, triangle, sdf };

class collision_model_t
{
  public:
    using volume_type = Eigen::AlignedBox3d; // In the future, enable different englobing volumes

    collision_model_t() = default;

    collision_model_t(collision_model_t const& other) = default;
    collision_model_t(collision_model_t&& other)      = default;

    collision_model_t& operator=(collision_model_t const& other) = default;
    collision_model_t& operator=(collision_model_t&& other) = default;

    virtual model_type_t model_type() const                                    = 0;
    virtual primitive_type_t primitive_type() const                            = 0;
    virtual void collide(collision_model_t& other, contact_handler_t& handler) = 0;
    virtual void update(simulation_t const& simulation)                        = 0;

    volume_type const& volume() const;
    volume_type& volume();

    index_type id() const;
    index_type& id();

    virtual ~collision_model_t() = default;

  private:
    volume_type englobing_volume_;
    index_type id_;
};

} // namespace collision
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_COLLISION_MODEL_H