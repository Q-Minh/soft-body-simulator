#ifndef SBS_PHYSICS_MODELS_LINEAR_TETRAHEDRAL_FEM_BODY_H
#define SBS_PHYSICS_MODELS_LINEAR_TETRAHEDRAL_FEM_BODY_H

#include "sbs/common/geometry.h"
#include "sbs/physics/body/body.h"
#include "sbs/physics/collision/bvh_model.h"
#include "sbs/physics/mechanics/linear_tetrahedral_fem_model.h"
#include "sbs/physics/visual/tetrahedral_fem_embedded_surface.h"

namespace sbs {
namespace physics {
namespace body {

class linear_tetrahedral_fem_body_t : public body_t
{
  public:
    using base_type = body_t;

    linear_tetrahedral_fem_body_t(
        xpbd::simulation_t& simulation,
        index_type const id,
        common::geometry_t const& geometry);

    virtual visual_model_type const& visual_model() const override;
    virtual collision_model_type const& collision_model() const override;
    virtual visual_model_type& visual_model() override;
    virtual collision_model_type& collision_model() override;
    virtual void update_visual_model() override;
    virtual void update_collision_model() override;
    virtual void update_physical_model() override;

  private:
    mechanics::linear_tetrahedral_fem_model_t mechanical_model_;
    visual::tetrahedral_fem_embedded_surface visual_model_;
    collision::point_bvh_model_t collision_model_;
};

} // namespace body
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MODELS_LINEAR_TETRAHEDRAL_FEM_BODY_H