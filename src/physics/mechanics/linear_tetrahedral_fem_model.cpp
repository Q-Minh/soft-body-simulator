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
        Eigen::Vector3d& xi       = this->dof(i);
        xi                        = Xi;
    }
}

typename linear_tetrahedral_fem_model_t::interpolation_function_type
linear_tetrahedral_fem_model_t::interpolation_field_at(Eigen::Vector3d const& X)
{
    index_type const e = domain().in_tetrahedron(X);
    assert(e != std::numeric_limits<index_type>::max());
    return interpolation_function_type(e, &(this->cells()), &(this->dofs()));
}

namespace differentiable {

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

} // namespace differentiable
} // namespace mechanics
} // namespace physics
} // namespace sbs