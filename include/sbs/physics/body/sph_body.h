#ifndef SBS_PHYSICS_BODY_SPH_BODY_H
#define SBS_PHYSICS_BODY_SPH_BODY_H

#include "body.h"
#include "sbs/physics/collision/bvh_model.h"
#include "sbs/physics/mechanics/sph_meshless_model.h"
#include "sbs/physics/visual/meshless_embedded_surface.h"

namespace sbs {
namespace physics {
namespace body {

template <class KernelFunctionType>
class sph_body_t : body_t
{
  public:
    using kernel_function_type  = KernelFunctionType;
    using mechanical_model_type = mechanics::sph_meshless_model_t<kernel_function_type>;

  private:
    mechanical_model_type mechanical_model_;
    visual::meshless_embedded_surface_t<mechanical_model_type> visual_model_;
    collision::point_bvh_model_t collision_model_;
};

} // namespace body
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_BODY_SPH_BODY_H