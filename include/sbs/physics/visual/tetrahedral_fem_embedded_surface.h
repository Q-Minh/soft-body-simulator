#ifndef SBS_PHYSICS_VISUAL_TETRAHEDRAL_FEM_EMBEDDED_SURFACE_H
#define SBS_PHYSICS_VISUAL_TETRAHEDRAL_FEM_EMBEDDED_SURFACE_H

#include "interpolated_embedded_surface.h"
#include "sbs/math/interpolation.h"
#include "sbs/physics/mechanics/linear_tetrahedral_fem_model.h"

#include <Eigen/Core>
#include <vector>

namespace sbs {
namespace physics {
namespace visual {

class tetrahedral_fem_embedded_surface
    : public interpolated_embedded_surface_t<math::differentiable::interpolation_op_t<
          math::differentiable::polynomial_hat_basis_function_t<1>>>
{
  public:
    using basis_function_type   = math::differentiable::polynomial_hat_basis_function_t<1>;
    using interpolation_op_type = math::differentiable::interpolation_op_t<basis_function_type>;
    using base_type             = interpolated_embedded_surface_t<interpolation_op_type>;

    tetrahedral_fem_embedded_surface() = default;

    tetrahedral_fem_embedded_surface(
        std::vector<Eigen::Vector3d> const& points,
        std::vector<index_type> const& indices,
        mechanics::linear_tetrahedral_fem_model_t const* mechanical_model);

    // Accessors
    mechanics::linear_tetrahedral_fem_model_t const* mechanical_model() const;
    index_type cell_containing_vertex(index_type vi) const;

    // Mutators

    /**
     * @brief
     * Computes new vertex positions based on the mechanical model's
     * current interpolation field.
     */
    void update();

  private:
    mechanics::linear_tetrahedral_fem_model_t const* mechanical_model_;
    std::vector<index_type> cell_containing_vertex_;
};

} // namespace visual
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_VISUAL_TETRAHEDRAL_FEM_EMBEDDED_SURFACE_H