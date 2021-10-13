#ifndef SBS_PHYSICS_ENVIRONMENT_BODY_H
#define SBS_PHYSICS_ENVIRONMENT_BODY_H

#include <sbs/common/mesh.h>
#include <sbs/physics/body.h>
#include <sbs/physics/collision/sdf_model.h>

namespace sbs {

namespace common {
struct geometry_t;
} // namespace common

namespace physics {

class environment_body_t : public body_t
{
  public:
    environment_body_t(
        simulation_t& simulation,
        index_type id,
        common::geometry_t const& geometry,
        std::array<unsigned int, 3u> const& resolution = {10, 10, 10});

    virtual visual_model_type const& visual_model() const override;
    virtual collision_model_type const& collision_model() const override;
    virtual visual_model_type& visual_model() override;
    virtual collision_model_type& collision_model() override;
    virtual void update_visual_model() override;
    virtual void update_collision_model() override;
    virtual void update_physical_model() override;
    virtual void transform(Eigen::Affine3d const& affine) override;

  private:
    common::static_mesh_t visual_model_;
    collision::sdf_model_t collision_model_;
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_ENVIRONMENT_BODY_H