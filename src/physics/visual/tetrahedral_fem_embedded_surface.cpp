#include "sbs/physics/visual/tetrahedral_fem_embedded_surface.h"

namespace sbs {
namespace physics {
namespace visual {

sbs::physics::visual::tetrahedral_fem_embedded_surface::tetrahedral_fem_embedded_surface(
    std::vector<Eigen::Vector3d> const& points,
    std::vector<index_type> const& indices,
    mechanics::linear_tetrahedral_fem_model_t const& mechanical_model)
    : base_type(points, indices), mechanical_model_(mechanical_model)
{
    using dof_type  = typename mechanics::linear_tetrahedral_fem_model_t::dof_type;
    using cell_type = typename mechanics::linear_tetrahedral_fem_model_t::cell_type;

    std::vector<interpolation_op_type> interpolation_operators{};
    interpolation_operators.reserve(this->vertex_count());

    geometry::tetrahedral_domain_t const& domain = mechanical_model_.domain();
    for (auto const& X : this->reference_positions())
    {
        index_type const ti   = domain.in_tetrahedron(X);
        cell_type const& cell = mechanical_model_.cell(ti);

        std::vector<basis_function_type> phis{};
        phis.reserve(cell.node_count());

        std::vector<dof_type> xis{};
        xis.reserve(cell.node_count());

        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            index_type const i              = cell.node(r);
            autodiff::Vector3dual const& xi = mechanical_model_.point(i);
            basis_function_type const& phi  = cell.phi(r);

            xis.push_back(xi);
            phis.push_back(phi);
        }

        interpolation_op_type const interp_op(xis, phis);
        interpolation_operators.push_back(interp_op);
    }

    this->use_interpolation_operators(interpolation_operators);
}

mechanics::linear_tetrahedral_fem_model_t const&
tetrahedral_fem_embedded_surface::mechanical_model() const
{
    return mechanical_model_;
}

void tetrahedral_fem_embedded_surface::update()
{
    auto const& Xs = this->reference_positions();

    for (auto i = 0u; i < Xs.size(); ++i)
    {
        auto const& Xi          = Xs[i];
        auto const& interpolate = this->interpolation_operator(i);
        auto const xi           = interpolate(Xi);
        this->position(i)       = xi.cast<scalar_type>();
    }
}

} // namespace visual
} // namespace physics
} // namespace sbs
