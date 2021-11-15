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

    efg_tetrahedral_meshless_model_t() = default;
    efg_tetrahedral_meshless_model_t(
        geometry::tetrahedral_domain_t const& domain,
        geometry::grid_t const& grid,
        scalar_type support);

    // Accessors
    geometry::tetrahedral_domain_t const& domain() const { return domain_; }
    geometry::grid_t const& grid() const { return grid_; }
    std::size_t integration_point_count() const { return efg_integration_points_.size(); }

    Eigen::Vector3d const& integration_point(index_type const i) const
    {
        return efg_integration_points_[i];
    }
    interpolation_function_type const&
    interpolation_field_from_integration_point(index_type const i) const
    {
        return interpolation_fields_[i];
    }
    std::vector<index_type> const& neighbours_of_integration_point(index_type const i) const
    {
        return neighbours_of_integration_point_[i];
    }

    interpolation_function_type interpolation_field_at(Eigen::Vector3d const& X) const;

    // Mutators

    /**
     * @brief
     * Adds an integration point Xi if it is in the tetrahedral domain.
     * Also adds its associated interpolation field.
     * @param Xi
     * @return true if the integration point and its associated interpolation field were added.
     * false otherwise
     */
    bool add_integration_point(Eigen::Vector3d const& Xi);

  private:
    geometry::tetrahedral_domain_t domain_; ///< The integration domain
    geometry::grid_t grid_; ///< The sampling grid used to sample meshless nodes and basis functions
                            ///< in the domain

    std::vector<Eigen::Vector3d>
        efg_integration_points_; ///< The integration points over the tetrahedral domain
    std::vector<interpolation_function_type>
        interpolation_fields_; ///< Interpolation fields around the integration points
    std::vector<std::vector<index_type>>
        neighbours_of_integration_point_; ///< Precomputed neighbourhoods
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
    this->set_support_radius(h);
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
}

template <class KernelFunctionType, unsigned int Order>
inline typename efg_tetrahedral_meshless_model_t<KernelFunctionType, Order>::
    interpolation_function_type
    efg_tetrahedral_meshless_model_t<KernelFunctionType, Order>::interpolation_field_at(
        Eigen::Vector3d const& X) const
{
    std::vector<index_type> const nodes = this->in_support_of_nodes(X);
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
    interpolation_op.cache_grad_phis(X);

    return interpolation_op;
}

template <class KernelFunctionType, unsigned int Order>
inline bool efg_tetrahedral_meshless_model_t<KernelFunctionType, Order>::add_integration_point(
    Eigen::Vector3d const& Xi)
{
    if (!domain_.contains(Xi))
        return false;

    Eigen::Vector3d const& integration_point = Xi;
    this->efg_integration_points_.push_back(integration_point);

    // Precompute the interpolation field around the integration point
    interpolation_function_type const interpolation_op = interpolation_field_at(integration_point);
    this->interpolation_fields_.push_back(interpolation_op);

    // Precompute neighbourhoods also
    std::vector<index_type> const nodes = this->in_support_of_nodes(integration_point);
    this->neighbours_of_integration_point_.push_back(nodes);

    return true;
}

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_EFG_TETRAHEDRAL_MESHLESS_MODEL_H