#ifndef SBS_PHYSICS_MECHANICS_EFG_TETRAHEDRAL_MESHLESS_MODEL_H
#define SBS_PHYSICS_MECHANICS_EFG_TETRAHEDRAL_MESHLESS_MODEL_H

#include "sbs/geometry/grid.h"
#include "sbs/geometry/tetrahedral_domain.h"
#include "sbs/math/meshless_model.h"
#include "sbs/math/mls.h"

#include <algorithm>

namespace sbs {
namespace physics {
namespace mechanics {

template <class KernelFunctionType, unsigned int Order>
class efg_tetrahedral_meshless_model_t : public math::meshless_model_t<
                                             autodiff::Vector3dual,
                                             math::mls_basis_function_t<KernelFunctionType, Order>>
{
  public:
    using kernel_function_type        = KernelFunctionType;
    using basis_function_type         = math::mls_basis_function_t<kernel_function_type, Order>;
    using interpolation_function_type = math::mls_interpolation_op_t<kernel_function_type, Order>;
    using self_type = efg_tetrahedral_meshless_model_t<kernel_function_type, Order>;

    efg_tetrahedral_meshless_model_t(
        geometry::tetrahedral_domain_t const& domain,
        geometry::grid_t const& grid,
        scalar_type support);

    geometry::tetrahedral_domain_t const& domain() const { return domain_; }
    geometry::grid_t const& grid() const { return grid_; }

    interpolation_function_type interpolation_field_at(Eigen::Vector3d const& X) const;

  private:
    geometry::tetrahedral_domain_t domain_; ///< The integration domain
    geometry::grid_t grid_; ///< The sampling grid used to sample meshless nodes and basis functions
                            ///< in the domain

    std::vector<Eigen::Vector3d>
        efg_integration_points_; ///< The integration points over the tetrahedral domain
    std::vector<interpolation_function_type>
        interpolation_fields_; ///< Interpolation fields around the integration points
};

template <class KernelFunctionType, unsigned int Order>
inline efg_tetrahedral_meshless_model_t<KernelFunctionType, Order>::
    efg_tetrahedral_meshless_model_t(
        geometry::tetrahedral_domain_t const& domain,
        geometry::grid_t const& grid,
        scalar_type support)
    : domain_(domain), grid_(grid), efg_integration_points_(), interpolation_fields_()
{
    // Sample meshless nodes
    grid_.walk([this](Eigen::Vector3d const& X) {
        bool const is_grid_cell_center_in_domain = domain_.contains(X);
        if (!is_grid_cell_center_in_domain)
            return;

        this->add_point(X);
        this->add_dof(X);
    });

    Eigen::Vector3d const& dX = grid_.delta();
    scalar_type const length  = dX.norm();
    scalar_type const h       = support * length;
    this->initialize_in_support_query(
        length * 1e-2); // Add a 1% tolerance to the spatial search data structure based on the
                        // grid cells' dimensions

    // Create node shape functions
    for (auto i = 0u; i < this->point_count(); ++i)
    {
        autodiff::Vector3dual const& Xi = this->point(i);
        std::vector<index_type> nodes   = this->in_support_of_nodes(Xi.cast<scalar_type>());
        kernel_function_type Wi(Xi, h);
        std::vector<autodiff::Vector3dual> Xjs{};
        std::vector<kernel_function_type> Wjs{};
        Xjs.reserve(nodes.size());
        Wjs.reserve(nodes.size());

        std::transform(
            nodes.begin(),
            nodes.end(),
            std::back_inserter(Xjs),
            [this](index_type const ni) {
                autodiff::Vector3dual const Xj = this->point(ni);
                return Xj;
            });
        std::transform(
            nodes.begin(),
            nodes.end(),
            std::back_inserter(Wjs),
            [this, h](index_type const ni) {
                autodiff::Vector3dual const Xj = this->point(ni);
                return kernel_function_type(Xj, h);
            });

        basis_function_type phi(Xi, Wi, Xjs, Wjs);
        this->add_basis_function(phi);
    }

    auto const& topology         = domain.topology();
    auto const tetrahedron_count = topology.tetrahedron_count();
    for (index_type ti = 0u; ti < tetrahedron_count; ++ti)
    {
        topology::tetrahedron_t const& t = topology.tetrahedron(ti);
        Eigen::Vector3d const& p1        = domain.position(t.v1());
        Eigen::Vector3d const& p2        = domain.position(t.v2());
        Eigen::Vector3d const& p3        = domain.position(t.v3());
        Eigen::Vector3d const& p4        = domain.position(t.v4());

        // Insert integration point at the barycenter of the integration domain's tetrahedra
        Eigen::Vector3d const integration_point = 0.25 * (p1 + p2 + p3 + p4);
        this->efg_integration_points_.push_back(integration_point);

        // Precompute the interpolation field around the integration point
        interpolation_function_type const interpolation_op =
            interpolation_field_at(integration_point);
        this->interpolation_fields_.push_back(interpolation_op);
    }
}

template <class KernelFunctionType, unsigned int Order>
inline typename efg_tetrahedral_meshless_model_t<KernelFunctionType, Order>::
    interpolation_function_type
    efg_tetrahedral_meshless_model_t<KernelFunctionType, Order>::interpolation_field_at(
        Eigen::Vector3d const& X) const
{
    std::vector<index_type> nodes = this->in_support_of_nodes(X);
    std::vector<autodiff::Vector3dual> xjs{};
    std::vector<basis_function_type> phijs{};
    xjs.reserve(nodes.size());
    phijs.reserve(nodes.size());

    std::transform(
        nodes.begin(),
        nodes.end(),
        std::back_inserter(xjs),
        [this](index_type const nj) {
            autodiff::Vector3dual const& xj = this->dof(nj);
            return xj;
        });
    std::transform(
        nodes.begin(),
        nodes.end(),
        std::back_inserter(phijs),
        [this](index_type const nj) {
            basis_function_type const& phij = this->phi(nj);
            return phij;
        });

    interpolation_function_type interpolation_op(xjs, phijs);
    interpolation_op.cache_grad_phis(integration_point);

    return interpolation_op;
}

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_EFG_TETRAHEDRAL_MESHLESS_MODEL_H