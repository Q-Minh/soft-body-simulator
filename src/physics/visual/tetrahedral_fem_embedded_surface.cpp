#include "sbs/physics/visual/tetrahedral_fem_embedded_surface.h"

namespace sbs {
namespace physics {
namespace visual {

tetrahedral_fem_embedded_surface::tetrahedral_fem_embedded_surface(
    std::vector<Eigen::Vector3d> const& points,
    std::vector<index_type> const& indices,
    mechanics::linear_tetrahedral_fem_model_t const* mechanical_model)
    : base_type(points, indices), mechanical_model_(mechanical_model)
{
    using dof_type  = typename mechanics::linear_tetrahedral_fem_model_t::dof_type;
    using cell_type = typename mechanics::linear_tetrahedral_fem_model_t::cell_type;

    std::vector<interpolation_op_type> interpolation_ops{};
    interpolation_ops.reserve(this->vertex_count());

    auto const& Xs = this->reference_positions();
    cell_containing_vertex_.resize(Xs.size());

    geometry::tetrahedral_domain_t const& domain = mechanical_model_->domain();
    for (auto vi = 0u; vi < Xs.size(); ++vi)
    {
        auto const& X = Xs[vi];

        index_type const ti         = domain.in_tetrahedron(X);
        cell_containing_vertex_[vi] = ti;
        cell_type const& cell       = mechanical_model_->cell(ti);

        std::vector<basis_function_type> phis{};
        phis.reserve(cell.node_count());

        std::vector<dof_type> xis{};
        xis.reserve(cell.node_count());

        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            index_type const i              = cell.node(r);
            autodiff::Vector3dual const& xi = mechanical_model_->point(i);
            basis_function_type const& phi  = cell.phi(r);

            xis.push_back(xi);
            phis.push_back(phi);
        }

        interpolation_op_type const interp_op(xis, phis);
        interpolation_ops.push_back(interp_op);
    }

    this->use_interpolation_operators(interpolation_ops);
    this->update();
}

mechanics::linear_tetrahedral_fem_model_t const*
tetrahedral_fem_embedded_surface::mechanical_model() const
{
    return mechanical_model_;
}

index_type tetrahedral_fem_embedded_surface::cell_containing_vertex(index_type vi) const
{
    return cell_containing_vertex_[vi];
}

void tetrahedral_fem_embedded_surface::update()
{
    auto const& Xs = this->reference_positions();

    for (auto i = 0u; i < Xs.size(); ++i)
    {
        auto const& Xi        = Xs[i];
        auto& interpolate     = this->interpolation_operator(i);
        index_type const ti   = cell_containing_vertex_[i];
        auto const& cell      = mechanical_model_->cell(ti);
        auto const node_count = cell.node_count();
        assert(interpolate.uis.size() == node_count && interpolate.phis.size() == node_count);
        // Update coefficients of the interpolation operation
        // using the mechanical model's interpolation field coefficients
        for (auto r = 0u; r < node_count; ++r)
        {
            auto const global_node_idx = cell.node(r);
            interpolate.uis[r]         = mechanical_model_->dof(global_node_idx);
        }
        auto const xi     = interpolate(Xi);
        this->position(i) = xi.cast<scalar_type>();
    }

    this->compute_normals();
}

} // namespace visual
} // namespace physics
} // namespace sbs
