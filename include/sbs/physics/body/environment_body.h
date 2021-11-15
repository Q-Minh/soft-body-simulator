#ifndef SBS_PHYSICS_ENVIRONMENT_BODY_H
#define SBS_PHYSICS_ENVIRONMENT_BODY_H

#include "sbs/common/mesh.h"
#include "body.h"
#include "sbs/physics/collision/sdf_model.h"

namespace sbs {

namespace common {
struct geometry_t;
} // namespace common

namespace physics {
namespace body {

class environment_body_t : public body_t
{
  public:
    environment_body_t(
        physics::xpbd::simulation_t& simulation,
        index_type id,
        common::geometry_t const& geometry,
        Eigen::AlignedBox3d const& domain,
        std::array<unsigned int, 3u> const& resolution = {10, 10, 10});

    environment_body_t(
        xpbd::simulation_t& simulation,
        index_type id,
        common::geometry_t const& geometry,
        collision::sdf_model_t const& sdf_model);

    virtual visual_model_type const& visual_model() const override;
    virtual collision_model_type const& collision_model() const override;
    virtual visual_model_type& visual_model() override;
    virtual collision_model_type& collision_model() override;
    virtual void update_visual_model() override;
    virtual void update_collision_model() override;
    virtual void update_physical_model() override;
    void transform(Eigen::Affine3d const& affine);

    collision::sdf_model_t const& sdf() const;

  private:
    common::static_mesh_t visual_model_;
    collision::sdf_model_t collision_model_;
};

} // namespace body
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_ENVIRONMENT_BODY_H