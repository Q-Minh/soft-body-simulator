#include "sbs/physics/mechanics/linear_tetrahedral_fem_model.h"

namespace sbs {
namespace physics {
namespace mechanics {

linear_tetrahedral_fem_model_t::linear_tetrahedral_fem_model_t(
    geometry::tetrahedral_domain_t const& domain)
    : base_type(domain)
{
    auto const num_points = this->point_count();
    for (auto i = 0u; i < num_points; ++i)
    {
        Eigen::Vector3d const& Xi = this->point(i);
        autodiff::Vector3dual& xi = this->dof(i);
        xi                        = autodiff::Vector3dual(Xi);
    }
}

} // namespace mechanics
} // namespace physics
} // namespace sbs