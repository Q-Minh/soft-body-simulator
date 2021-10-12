#ifndef SBS_PHYSICS_ENVIRONMENT_BODY_H
#define SBS_PHYSICS_ENVIRONMENT_BODY_H

#include <sbs/common/mesh.h>
#include <sbs/physics/body.h>
#include <sbs/physics/collision/sdf_model.h>

namespace sbs {
namespace physics {

class environment_body_t : public body_t
{
  public:
    environment_body_t(
        common::geometry_t const& geometry,
        std::array<unsigned int, 3u> const& resolution = {10, 10, 10});

    virtual visual_model_type const& visual_model() const override;
    virtual collision_model_type const& collision_model() const override;
    virtual void update_visual_model(simulation_t const& simulation) override;
    virtual void update_collision_model(simulation_t const& simulation) override;
    virtual void update_physical_model(simulation_t const& simulation) override;

  private:
    common::static_mesh_t visual_model_;
    collision::sdf_model_t collision_model_;
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_ENVIRONMENT_BODY_H