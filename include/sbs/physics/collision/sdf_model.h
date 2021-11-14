#ifndef SBS_PHYSICS_COLLISION_SDF_MODEL_H
#define SBS_PHYSICS_COLLISION_SDF_MODEL_H

#include "sbs/aliases.h"
#include "sbs/physics/collision/collision_model.h"

#include <Discregrid/cubic_lagrange_discrete_grid.hpp>
#include <utility>

namespace sbs {
namespace physics {
namespace collision {

class contact_handler_t;

class sdf_model_t : public collision_model_t
{
  public:
    using analytic_sdf_type =
        std::function<std::pair<scalar_type, Eigen::Vector3d>(Eigen::Vector3d const&)>;

    sdf_model_t(Eigen::AlignedBox3d const& domain, std::array<unsigned int, 3u> const& resolution);
    sdf_model_t(Discregrid::CubicLagrangeDiscreteGrid const& sdf);
    sdf_model_t(analytic_sdf_type const& analytic_sdf, Eigen::AlignedBox3d const& volume);

    sdf_model_t(sdf_model_t const& other) = default;
    sdf_model_t(sdf_model_t&& other)      = default;

    sdf_model_t& operator=(sdf_model_t const& other) = default;
    sdf_model_t& operator=(sdf_model_t&& other) noexcept = default;

    virtual model_type_t model_type() const override;
    virtual primitive_type_t primitive_type() const override;
    virtual void collide(collision_model_t& other, contact_handler_t& handler) override;
    virtual void update(xpbd::simulation_t const& simulation) override;

    static sdf_model_t
    from_plane(Eigen::Hyperplane<scalar_type, 3> const& plane, Eigen::AlignedBox3d const& volume);

    std::pair<scalar_type, Eigen::Vector3d> evaluate(Eigen::Vector3d const& p) const;

  private:
    Discregrid::CubicLagrangeDiscreteGrid sdf_;
    analytic_sdf_type analytic_sdf_;
};

} // namespace collision
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_SDF_MODEL_H