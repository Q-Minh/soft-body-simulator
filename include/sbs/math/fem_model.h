#ifndef SBS_MATH_FEM_MODEL_H
#define SBS_MATH_FEM_MODEL_H

#include "basis_functions.h"
#include "cell.h"
#include "element.h"
#include "sbs/geometry/tetrahedral_domain.h"

#include <Eigen/Core>
#include <vector>

namespace sbs {
namespace math {

/**
 * @brief
 * A Finite Element Method data structure containing
 * all necessary information for solving PDEs using FEM
 * on arbitrary elements with arbitrary basis functions
 * of the given order. The unknowns are stored in the dofs
 * of this model. Assumes a Galerkin method where trial and
 * test functions are in the same space.
 * @tparam DofType Type of the unknowns to solve for
 * @tparam CellType Type of the cell containing the nodes
 * @tparam ElementType Type of the elements of the FEM domain
 */
template <class DofType, class ElementType, class CellType>
class fem_model_t
{
  public:
    using dof_type     = DofType;
    using cell_type    = CellType;
    using element_type = ElementType;
    using size_type    = std::size_t;

    fem_model_t() = default;

    // Accessors
    dof_type const& dof(index_type i) const { return dofs_[i]; }
    autodiff::Vector3dual const& point(index_type i) const { return points_[i]; }
    element_type const& element(index_type e) const { return elements_[e]; }
    cell_type const& cell(index_type e) const { return cells_[e]; }

    dof_type& dof(index_type i) { return dofs_[i]; }

    size_type dof_count() const { return dofs_.size(); }
    size_type point_count() const { return points_.size(); }
    size_type element_count() const { return elements_.size(); }
    size_type cell_count() const { return cells_.size(); }

    // Modifiers
    void add_dof(dof_type const& dof) { dofs_.push_back(dof); }
    void add_point(autodiff::Vector3dual const& point) { points_.push_back(point); }
    void add_element(element_type const& element) { elements_.push_back(element); }
    void add_cell(cell_type const& cell) { cells_.push_back(cell); }

    void clear_dofs() { dofs_.clear(); }
    void clear_points() { points_.clear(); }
    void clear_elements() { elements_.clear(); }
    void clear_cells() { cells_.clear(); }
    void clear()
    {
        clear_dofs();
        clear_points();
        clear_elements();
        clear_cells();
    }

  private:
    std::vector<dof_type> dofs_;
    std::vector<autodiff::Vector3dual> points_;
    std::vector<element_type> elements_;
    std::vector<cell_type> cells_;
};

/**
 * @brief
 * A Finite Element Method data structure containing
 * all necessary information for solving PDEs using FEM
 * on tetrahedral elements with polynomial hat basis functions
 * of the given order. The unknowns are stored in the dofs
 * of this model.
 * @tparam DofType Type of the unknown coefficients to solve for
 */
template <class DofType, unsigned int Order>
class tetrahedral_fem_model_t : public fem_model_t<
                                    DofType,
                                    tetrahedral_element_t,
                                    cell_t<
                                        num_nodes_for_tetrahedron_cell_of_order_t<Order>::value,
                                        polynomial_hat_basis_function_t<Order>>>
{
  public:
    using base_type = fem_model_t<
        DofType,
        tetrahedral_element_t,
        cell_t<
            num_nodes_for_tetrahedron_cell_of_order_t<Order>::value,
            polynomial_hat_basis_function_t<Order>>>;

    using dof_type = DofType;

    using cell_type = cell_t<
        num_nodes_for_tetrahedron_cell_of_order_t<Order>::value,
        polynomial_hat_basis_function_t<Order>>;

    using element_type        = tetrahedral_element_t;
    using basis_function_type = polynomial_hat_basis_function_t<Order>;
    using self_type           = tetrahedral_fem_model_t<dof_type, Order>;

    tetrahedral_fem_model_t() = default;
    tetrahedral_fem_model_t(geometry::tetrahedral_domain_t const& domain);

    // Accessors
    geometry::tetrahedral_domain_t const& domain() const { return domain_; }

    // Mutators
    void build_model();

  private:
    geometry::tetrahedral_domain_t domain_;
};

template <class DofType, unsigned int Order>
tetrahedral_fem_model_t<DofType, Order>::tetrahedral_fem_model_t(
    geometry::tetrahedral_domain_t const& domain)
    : domain_(domain)
{
    build_model();
}

template <class DofType, unsigned int Order>
void tetrahedral_fem_model_t<DofType, Order>::build_model()
{
    this->clear();

    auto const& topology     = domain_.topology();
    auto const element_count = topology.tetrahedron_count();
    auto const vertex_count  = topology.vertex_count();

    // Create elements, cells, points and dofs by traversing the domain's tetrahedra
    for (auto e = 0u; e < element_count; ++e)
    {
        topology::tetrahedron_t const& tetrahedron = topology.tetrahedron(e);
        autodiff::Vector3dual const X1             = domain_.position(tetrahedron.v1());
        autodiff::Vector3dual const X2             = domain_.position(tetrahedron.v2());
        autodiff::Vector3dual const X3             = domain_.position(tetrahedron.v3());
        autodiff::Vector3dual const X4             = domain_.position(tetrahedron.v4());

        element_type element(X1, X2, X3, X4);
        unsigned int constexpr node_count = cell_type::node_count_value;

        std::vector<autodiff::Vector3dual> Xis{};
        Xis.reserve(node_count);

        // Sample nodes uniformly in the tetrahedral element
        unsigned int const num_segments_on_edge = Order;
        scalar_type deltaX = 1. / static_cast<scalar_type>(num_segments_on_edge);
        unsigned int constexpr num_samples_on_edge = Order + 1u;
        for (auto k = 0u; k < num_samples_on_edge; ++k)
        {
            auto const j_samples = num_samples_on_edge - k;
            for (auto j = 0u; j < j_samples; ++j)
            {
                auto const i_samples = num_samples_on_edge - k - j;
                for (auto i = 0u; i < i_samples; ++i)
                {
                    auto const dX = i * deltaX;
                    auto const dY = j * deltaX;
                    auto const dZ = k * deltaX;

                    autodiff::Vector3dual const Xi = X1 + autodiff::Vector3dual(dX, dY, dZ);
                    Xis.push_back(Xi);
                }
            }
        }

        // Build the polynomial hat basis functions for each sampled node
        Eigen::Matrix<autodiff::dual, node_count, node_count> P{};
        for (auto i = 0u; i < node_count; ++i)
        {
            P.row(i) = polynomial3d<Order>(Xis[i]);
        }

        using PolynomialMatrixType      = Eigen::Matrix<autodiff::dual, node_count, node_count>;
        PolynomialMatrixType const Pinv = P.inverse(); // LU-decomposition based inverse

        auto const node_index_offset = this->point_count();
        cell_type cell{};
        for (auto r = 0u; r < node_count; ++r)
        {
            basis_function_type const phi(Pinv.col(r));
            autodiff::Vector3dual const& Xr = Xis[r];
            index_type const i              = static_cast<index_type>(node_index_offset + r);

            cell.set_node(r, i);
            cell.set_phi(r, phi);
        }

        // build fem model
        this->add_element(element);
        this->add_cell(cell);
        for (auto r = 0u; r < node_count; ++r)
        {
            this->add_point(Xis[r]);
            this->add_dof(dof_type{});
        }
    }
}

} // namespace math
} // namespace sbs

#endif // SBS_MATH_FEM_MODEL_H