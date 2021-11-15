#ifndef SBS_PHYSICS_MECHANICS_LINEAR_TETERAHEDRAL_FEM_MODEL_H
#define SBS_PHYSICS_MECHANICS_LINEAR_TETERAHEDRAL_FEM_MODEL_H

#include "sbs/geometry/tetrahedral_domain.h"
#include "sbs/math/fem_model.h"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace sbs {
namespace physics {
namespace mechanics {

/**
 * @brief
 * Solid elastic body modeled with FEM using linear polynomial basis function in tetrahedral
 * elements.
 */
class linear_tetrahedral_fem_model_t : public math::tetrahedral_fem_model_t<autodiff::Vector3dual, 1>
{
  public:
    using base_type = math::tetrahedral_fem_model_t<autodiff::Vector3dual, 1>;

    linear_tetrahedral_fem_model_t() = default;

    /**
     * @brief
     * Initializes the fem model with a tetrahedral mesh as the simulated domain.
     * Initializes the degrees of freedom, i.e. world space positions, associated
     * wtih each point in the interpolation field of the FEM model.
     * @param domain The (tetrahedral) domain to simulate.
     */
    linear_tetrahedral_fem_model_t(geometry::tetrahedral_domain_t const& domain);

    autodiff::Vector3dual const& x(index_type i) const { return this->dof(i); }
    autodiff::Vector3dual const& X(index_type i) const { return this->point(i); }
    autodiff::Vector3dual u(index_type i) const { return x(i) - X(i); }

  private:
};

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_LINEAR_TETERAHEDRAL_FEM_MODEL_H