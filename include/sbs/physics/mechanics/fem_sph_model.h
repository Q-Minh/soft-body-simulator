#ifndef SBS_PHYSICS_MECHANICS_FEM_SPH_MODEL_H
#define SBS_PHYSICS_MECHANICS_FEM_SPH_MODEL_H

#include "sbs/math/fem_model.h"
#include "sbs/math/fem_sph.h"
#include "sbs/math/meshless_model.h"

namespace sbs {
namespace physics {
namespace mechanics {

template <class KernelFunctionType>
class fem_sph_model_t : math::tetrahedral_fem_model_t<Eigen::Vector3d, 1u>
{
  public:
  private:
};

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_FEM_SPH_MODEL_H
