#ifndef SBS_PHYSICS_COLLISION_SDF_MODEL_H
#define SBS_PHYSICS_COLLISION_SDF_MODEL_H

#include <Discregrid/cubic_lagrange_discrete_grid.hpp>
#include <sbs/aliases.h>
#include <sbs/physics/collision/collision_model.h>
#include <utility>

namespace sbs {
namespace physics {
namespace collision {

class contact_handler_t;

class sdf_model_t : public collision_model_t
{
  public:
    sdf_model_t(Eigen::AlignedBox3d const& domain, std::array<unsigned int, 3u> const& resolution);
    sdf_model_t(Discregrid::CubicLagrangeDiscreteGrid const& sdf);

    sdf_model_t(sdf_model_t const& other) = default;
    sdf_model_t(sdf_model_t&& other)      = default;

    sdf_model_t& operator=(sdf_model_t const& other) = default;
    sdf_model_t& operator=(sdf_model_t&& other) = default;

    virtual model_type_t model_type() const override;
    virtual primitive_type_t primitive_type() const override;
    virtual void collide(collision_model_t& other, contact_handler_t& handler) override;
    virtual void update(simulation_t const& simulation) override;

    static sdf_model_t from_plane(
        Eigen::AlignedBox3d const& domain,
        std::array<unsigned int, 3u> const& resolution,
        Eigen::Hyperplane<scalar_type, 3> const& plane,
        scalar_type const sampling_depth_from_surface = scalar_type{1.});

    std::pair<scalar_type, Eigen::Vector3d> evaluate(Eigen::Vector3d const& p) const;

  private:
    Discregrid::CubicLagrangeDiscreteGrid sdf_;
};

} // namespace collision
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_SDF_MODEL_H