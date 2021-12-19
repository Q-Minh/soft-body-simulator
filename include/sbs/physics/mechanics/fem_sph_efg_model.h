#ifndef SBS_PHYSICS_MECHANICS_FEM_SPH_EFG_MODEL_H
#define SBS_PHYSICS_MECHANICS_FEM_SPH_EFG_MODEL_H

#include "sbs/geometry/grid.h"
#include "sbs/math/fem.h"
#include "sbs/math/fem_model.h"
#include "sbs/math/fem_sph.h"
#include "sbs/math/interpolation.h"
#include "sbs/math/meshless_model.h"

namespace sbs {
namespace physics {
namespace mechanics {

template <class KernelFunctionType>
class fem_sph_efg_model_t : public math::tetrahedral_fem_model_t<Eigen::Vector3d, 1u, false>
{
  public:
    using base_type = math::tetrahedral_fem_model_t<Eigen::Vector3d, 1u, false>;

    using dof_type             = Eigen::Vector3d;
    using point_type           = Eigen::Vector3d;
    using kernel_function_type = KernelFunctionType;

    using fem_basis_function_type = basis_function_type;
    using sph_basis_function_type =
        math::differentiable::sph_basis_function_t<kernel_function_type>;

    using fem_model_type = math::tetrahedral_fem_model_t<Eigen::Vector3d, 1u, false>;
    using sph_model_type = math::meshless_model_t<dof_type, point_type, sph_basis_function_type>;
    using meshless_model_type = sph_model_type;

    using fem_interpolation_function_type =
        math::fem_interpolation_t<typename base_type::cell_type>;

    using mixed_interpolation_function_type =
        math::fem_sph_efg_interpolation_t<typename base_type::cell_type, kernel_function_type>;

    using mixed_deformation_gradient_function_type = math::
        fem_sph_efg_deformation_gradient_op_t<typename base_type::cell_type, kernel_function_type>;

    using self_type = fem_sph_efg_model_t<KernelFunctionType>;

    fem_sph_efg_model_t() = default;
    fem_sph_efg_model_t(
        geometry::tetrahedral_domain_t const& domain,
        geometry::grid_t const& grid,
        scalar_type support);

    fem_sph_efg_model_t(self_type const& other);
    fem_sph_efg_model_t& operator=(self_type const& other);

    // Accessors
    sph_model_type const& meshless_model() const { return sph_model_; }
    sph_model_type& meshless_model() { return sph_model_; }

    scalar_type const& V(index_type j) const { return Vjs_[j]; }
    kernel_function_type const& W(index_type j) const { return Wjs_[j]; }
    Eigen::Vector3d const integration_point(index_type k) const { return Xks_[k]; }
    index_type const tetrahedron_around_integration_point(index_type k) const
    {
        return efg_point_tets_[k];
    }
    index_type const integration_point_in_tetrahedron(index_type t) const { return tet_to_efg_[t]; }

    fem_interpolation_function_type const& fem_interpolation_field_at(index_type e) const
    {
        return fem_interpolation_fields_[e];
    }
    mixed_interpolation_function_type const& mixed_interpolation_field_at(index_type k) const
    {
        return mixed_interpolation_fields_[k];
    }
    mixed_interpolation_function_type& mixed_interpolation_field_at(index_type k)
    {
        return mixed_interpolation_fields_[k];
    }
    mixed_interpolation_function_type mixed_interpolation_field_at(Eigen::Vector3d const& X);

    mixed_deformation_gradient_function_type const&
    mixed_deformation_gradient_function(index_type k) const
    {
        return mixed_deformation_gradient_functions_[k];
    }
    bool is_mixed_cell(index_type const e) const
    {
        return tet_to_efg_[e] != std::numeric_limits<index_type>::max();
    }

    std::vector<index_type> const& neighbours_of_meshless_node(index_type const j) const
    {
        return mixed_interpolation_fields_[j].js;
    }

    bool has_basis_function(index_type const i) const { return has_basis_function_[i]; }

    std::size_t num_active_nodes() const
    {
        return std::count(has_basis_function_.begin(), has_basis_function_.end(), true);
    }
    std::size_t num_sph_nodes() const { return sph_model_.dof_count(); }
    std::size_t total_dof_count() const { return num_active_nodes() + sph_model_.dof_count(); }
    std::size_t efg_integration_point_count() const { return Xks_.size(); }

  private:
    sph_model_type sph_model_; ///< SPH particles in the tetrahedral domain

    std::vector<scalar_type> Vjs_;          ///< Meshless shepard coefficients (nodal volume)
    std::vector<kernel_function_type> Wjs_; ///< Meshless nodal kernel functions

    std::vector<Eigen::Vector3d> Xks_; ///< EFG integration points
    std::vector<mixed_interpolation_function_type>
        mixed_interpolation_fields_; ///< Interpolation fields at EFG points
    std::vector<mixed_deformation_gradient_function_type>
        mixed_deformation_gradient_functions_; ///< Deformation gradient functions at EFG points

    std::vector<std::vector<index_type>> particles_in_tet_; ///< Particles in each tetrahedron
    std::vector<bool> has_basis_function_; ///< Marks active fem dofs vs inactive fem dofs

    std::vector<fem_interpolation_function_type>
        fem_interpolation_fields_; ///< Interpolations in interior tetrahedra

    std::vector<index_type> efg_point_tets_; ///< The englobing tetrahedron of each efg point
    std::vector<index_type>
        tet_to_efg_; ///< Gives index of EFG integration in the given tetrahedron
};

template <class KernelFunctionType>
inline fem_sph_efg_model_t<KernelFunctionType>::fem_sph_efg_model_t(
    geometry::tetrahedral_domain_t const& domain,
    geometry::grid_t const& grid,
    scalar_type support)
    : base_type(domain),
      sph_model_(),
      Vjs_(),
      Wjs_(),
      Xks_(),
      mixed_interpolation_fields_(),
      mixed_deformation_gradient_functions_(),
      has_basis_function_(),
      fem_interpolation_fields_(),
      efg_point_tets_(),
      tet_to_efg_()
{
    // Initialize fem interpolation fields. mixed cells should not use their fem interpolation
    // fields
    this->fem_interpolation_fields_.reserve(this->cell_count());
    for (auto e = 0u; e < this->cell_count(); ++e)
    {
        fem_interpolation_function_type const fem_interpolation(
            e,
            &(this->cells()),
            &(this->dofs()));
        this->fem_interpolation_fields_.push_back(fem_interpolation);
    }

    // Sample meshless nodes at the boundary tet vertices
    auto const& topology    = this->domain().topology();
    auto const vertex_count = topology.vertex_count();
    has_basis_function_.resize(vertex_count, true);
    auto const interior_vertex_indices = topology.interior_vertex_indices();
    for (index_type const vi : interior_vertex_indices)
    {
        this->dof(vi) = this->point(vi);
    }
    auto const boundary_vertex_indices = topology.boundary_vertex_indices();
    for (index_type const vi : boundary_vertex_indices)
    {
        Eigen::Vector3d const& Xi = this->domain().position(vi);
        sph_model_.add_point(Xi);
        sph_model_.add_dof(Xi);
        has_basis_function_[vi] = false;
    }

    // Initialize meshless shape function support with longest tetrahedral edge
    scalar_type max_edge_length = 0.;
    auto const edge_count       = topology.edge_count();
    for (auto i = 0u; i < edge_count; ++i)
    {
        index_type const ei           = static_cast<index_type>(i);
        topology::edge_t const& e     = topology.edge(ei);
        Eigen::Vector3d const& X1     = this->domain().position(e.v1());
        Eigen::Vector3d const& X2     = this->domain().position(e.v2());
        scalar_type const edge_length = (X2 - X1).norm();
        if (edge_length > max_edge_length)
        {
            max_edge_length = edge_length;
        }
    }
    scalar_type const h = support * max_edge_length;
    sph_model_.set_support_radius(h);
    sph_model_.initialize_in_support_query(max_edge_length * 1e-2);

    // Create nodal shape functions and neighbourhoods
    std::vector<scalar_type> Vjs{};
    Vjs.reserve(sph_model_.point_count());
    std::vector<kernel_function_type> Wjs{};
    Wjs.reserve(sph_model_.point_count());
    std::vector<std::vector<index_type>> Njks{};
    Njks.reserve(sph_model_.point_count());

    for (auto j = 0u; j < sph_model_.point_count(); ++j)
    {
        Eigen::Vector3d const& Xj = sph_model_.point(j);
        kernel_function_type Wj(Xj, h);
        Wjs.push_back(Wj);

        std::vector<index_type> nodes = sph_model_.in_support_of_nodes(Xj);
        Njks.push_back(nodes);
    }
    for (auto j = 0u; j < sph_model_.point_count(); ++j)
    {
        Eigen::Vector3d const& Xj          = sph_model_.point(j);
        std::vector<index_type> const& Nks = Njks[j];
        scalar_type const density =
            std::accumulate(Nks.begin(), Nks.end(), 0., [&](scalar_type sum, index_type const Nk) {
                kernel_function_type const& Wk = Wjs[Nk];
                auto const dualW               = Wk(Xj);
                scalar_type const W            = static_cast<scalar_type>(dualW);
                return sum + W;
            });
        scalar_type const Vj = 1. / density;
        Vjs.push_back(Vj);
        kernel_function_type const& Wj = Wjs[j];
        sph_basis_function_type const sph_basis_function(Wj, Vj);
        sph_model_.add_basis_function(sph_basis_function);
    }

    this->Vjs_ = Vjs;
    this->Wjs_ = Wjs;

    // Sample EFG points, interpolations and gradients
    auto const boundary_tetrahedron_indices = topology.boundary_tetrahedron_indices();
    Xks_.reserve(boundary_tetrahedron_indices.size());
    efg_point_tets_.reserve(boundary_tetrahedron_indices.size());
    tet_to_efg_.resize(topology.tetrahedron_count(), std::numeric_limits<index_type>::max());
    for (index_type const ti : boundary_tetrahedron_indices)
    {
        index_type const k = static_cast<index_type>(Xks_.size());

        topology::tetrahedron_t const& t = topology.tetrahedron(ti);
        Eigen::Vector3d const X1         = this->domain().position(t.v1());
        Eigen::Vector3d const X2         = this->domain().position(t.v2());
        Eigen::Vector3d const X3         = this->domain().position(t.v3());
        Eigen::Vector3d const X4         = this->domain().position(t.v4());

        Eigen::Vector3d const Xk = 0.25 * (X1 + X2 + X3 + X4);
        Xks_.push_back(Xk);
        efg_point_tets_.push_back(ti);
        tet_to_efg_[ti] = k;

        auto const& cell                 = this->cell(ti);
        std::vector<index_type> const js = sph_model_.in_support_of_nodes(Xk);
        mixed_interpolation_function_type mixed_interpolation(
            Xk,
            ti,
            &(this->cells()),
            &(this->points()),
            &(this->dofs()),
            &(this->has_basis_function_),
            js,
            &(sph_model_.points()),
            &(sph_model_.dofs()),
            &(this->Vjs_),
            &(this->Wjs_));

        this->mixed_interpolation_fields_.push_back(mixed_interpolation);

        mixed_deformation_gradient_function_type mixed_deformation_gradient(
            &(this->Xks_),
            k,
            this->mixed_interpolation_fields_[k]);
        this->mixed_deformation_gradient_functions_.push_back(mixed_deformation_gradient);
    }
}

template <class KernelFunctionType>
inline fem_sph_efg_model_t<KernelFunctionType>::fem_sph_efg_model_t(self_type const& other)
    : base_type(domain),
      sph_model_(other.sph_model_),
      Vjs_(other.Vjs_),
      Wjs_(other.Wjs_),
      Xks_(other.Xks_),
      mixed_interpolation_fields_(),
      mixed_deformation_gradient_functions_(),
      has_basis_function_(other.has_basis_function_),
      fem_interpolation_fields_(),
      efg_point_tets_(other.efg_point_tets_),
      tet_to_efg_(other.tet_to_efg_)
{
    mixed_interpolation_fields_.reserve(other.efg_integration_point_count());
    mixed_deformation_gradient_functions_.reserve(other.efg_integration_point_count());
    for (auto k = 0u; k < other.efg_integration_point_count(); ++k)
    {
        Eigen::Vector3d const& Xk = other.integration_point(k);
        mixed_interpolation_function_type const& other_interpolation =
            other.mixed_interpolation_fields_[k];
        index_type const ti = efg_point_tets_[k];
        mixed_interpolation_function_type interpolation(
            Xk,
            ti,
            &(this->cells()),
            &(this->points()),
            &(this->dofs()),
            &(this->has_basis_function_),
            other_interpolation.js,
            &(sph_model_.points()),
            &(sph_model_.dofs()),
            &(this->Vjs_),
            &(this->Wjs_));
        mixed_interpolation_fields_.push_back(interpolation);

        mixed_deformation_gradient_function_type deformation_gradient_function(
            &(this->Xks_),
            k,
            this->mixed_interpolation_fields_[k]);
        this->mixed_deformation_gradient_functions_.push_back(deformation_gradient_function);
    }

    this->fem_interpolation_fields_.reserve(this->cell_count());
    for (auto e = 0u; e < this->cell_count(); ++e)
    {
        fem_interpolation_function_type const fem_interpolation(
            e,
            &(this->cells()),
            &(this->dofs()));
        this->fem_interpolation_fields_.push_back(fem_interpolation);
    }
}

template <class KernelFunctionType>
inline fem_sph_efg_model_t<KernelFunctionType>&
fem_sph_efg_model_t<KernelFunctionType>::operator=(self_type const& other)
{
    base_type::operator=(other);

    sph_model_          = other.sph_model_;
    Vjs_                = other.Vjs_;
    Wjs_                = other.Wjs_;
    Xks_                = other.Xks_;
    has_basis_function_ = other.has_basis_function_;
    efg_point_tets_     = other.efg_point_tets_;
    tet_to_efg_         = other.tet_to_efg_;

    mixed_interpolation_fields_.reserve(other.efg_integration_point_count());
    mixed_deformation_gradient_functions_.reserve(other.efg_integration_point_count());
    for (auto k = 0u; k < other.efg_integration_point_count(); ++k)
    {
        Eigen::Vector3d const& Xk = other.integration_point(k);
        mixed_interpolation_function_type const& other_interpolation =
            other.mixed_interpolation_fields_[k];
        index_type const ti = efg_point_tets_[k];
        mixed_interpolation_function_type interpolation(
            Xk,
            ti,
            &(this->cells()),
            &(this->points()),
            &(this->dofs()),
            &(this->has_basis_function_),
            other_interpolation.js,
            &(sph_model_.points()),
            &(sph_model_.dofs()),
            &(this->Vjs_),
            &(this->Wjs_));
        mixed_interpolation_fields_.push_back(interpolation);

        mixed_deformation_gradient_function_type deformation_gradient_function(
            &(this->Xks_),
            k,
            this->mixed_interpolation_fields_[k]);
        this->mixed_deformation_gradient_functions_.push_back(deformation_gradient_function);
    }

    this->fem_interpolation_fields_.reserve(this->cell_count());
    for (auto e = 0u; e < this->cell_count(); ++e)
    {
        fem_interpolation_function_type const fem_interpolation(
            e,
            &(this->cells()),
            &(this->dofs()));
        this->fem_interpolation_fields_.push_back(fem_interpolation);
    }

    return *this;
}

template <class KernelFunctionType>
inline typename fem_sph_efg_model_t<KernelFunctionType>::mixed_interpolation_function_type
fem_sph_efg_model_t<KernelFunctionType>::mixed_interpolation_field_at(Eigen::Vector3d const& X)
{
    index_type const ti              = this->domain().in_tetrahedron(X);
    std::vector<index_type> const js = sph_model_.in_support_of_nodes(X);

    mixed_interpolation_function_type interpolation(
        X,
        ti,
        &(this->cells()),
        &(this->points()),
        &(this->dofs()),
        &(this->has_basis_function_),
        js,
        &(sph_model_.points()),
        &(sph_model_.dofs()),
        &(this->Vjs_),
        &(this->Wjs_));

    return interpolation;
}

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_FEM_SPH_EFG_MODEL_H
