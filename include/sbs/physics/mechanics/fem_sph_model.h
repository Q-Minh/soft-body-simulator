#ifndef SBS_PHYSICS_MECHANICS_FEM_SPH_MODEL_H
#define SBS_PHYSICS_MECHANICS_FEM_SPH_MODEL_H

#include "sbs/geometry/grid.h"
#include "sbs/math/fem.h"
#include "sbs/math/fem_model.h"
#include "sbs/math/fem_sph.h"
#include "sbs/math/interpolation.h"
#include "sbs/math/meshless_model.h"

#include <numeric>

namespace sbs {
namespace physics {
namespace mechanics {

template <class KernelFunctionType>
class fem_sph_model_t
    : public math::tetrahedral_fem_model_t<Eigen::Vector3d, 1u, false /*Don't use autodiff types*/>
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
        math::fem_sph_interpolation_t<typename base_type::cell_type, kernel_function_type>;

    using mixed_deformation_gradient_function_type = math::fem_sph_nodal_deformation_gradient_op_t<
        typename base_type::cell_type,
        kernel_function_type>;

    using self_type = fem_sph_model_t<KernelFunctionType>;

    fem_sph_model_t() = default;
    fem_sph_model_t(
        geometry::tetrahedral_domain_t const& domain,
        geometry::grid_t const& grid,
        scalar_type support);

    fem_sph_model_t(self_type const& other);
    fem_sph_model_t& operator=(self_type const& other);

    // Accessors
    sph_model_type const& meshless_model() const { return sph_model_; }
    sph_model_type& meshless_model() { return sph_model_; }
    geometry::grid_t const& grid() const { return grid_; }

    scalar_type const& V(index_type j) const { return Vjs_[j]; }
    kernel_function_type const& W(index_type j) const { return Wjs_[j]; }
    Eigen::Matrix3d const& F(index_type j) const { return Fjs_[j]; }
    Eigen::Matrix3d& F(index_type j) { return Fjs_[j]; }

    fem_interpolation_function_type const& fem_interpolation_field_at(index_type e) const
    {
        return fem_interpolation_fields_[e];
    }
    mixed_interpolation_function_type const& mixed_interpolation_field_at(index_type j) const
    {
        return mixed_interpolation_fields_[j];
    }
    mixed_interpolation_function_type& mixed_interpolation_field_at(index_type j)
    {
        return mixed_interpolation_fields_[j];
    }
    mixed_interpolation_function_type mixed_interpolation_field_at(Eigen::Vector3d const& X);

    mixed_deformation_gradient_function_type const&
    mixed_deformation_gradient_function(index_type j) const
    {
        return mixed_deformation_gradient_functions_[j];
    }

    std::vector<index_type> const& particles_in_tetrahedron(index_type const ti) const
    {
        return particles_in_tet_[ti];
    }
    index_type englobing_tetrahedron_of_particle(index_type const j) const
    {
        return particle_tets_[j];
    }
    bool is_mixed_cell(index_type const e) const { return !particles_in_tet_[e].empty(); }

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

  private:
    sph_model_type sph_model_; ///< SPH particles in the tetrahedral domain

    geometry::grid_t grid_; ///< The particle sampling grid

    std::vector<scalar_type> Vjs_;          ///< Meshless shepard coefficients (nodal volume)
    std::vector<kernel_function_type> Wjs_; ///< Meshless nodal kernel functions
    std::vector<Eigen::Matrix3d> Fjs_;      ///< Deformation gradients at each meshless point

    std::vector<mixed_interpolation_function_type>
        mixed_interpolation_fields_; ///< Interpolation fields at SPH points
    std::vector<mixed_deformation_gradient_function_type>
        mixed_deformation_gradient_functions_; ///< SPH nodal deformation gradient functions

    std::vector<std::vector<index_type>> particles_in_tet_; ///< Particles in each tetrahedron
    std::vector<index_type> particle_tets_; ///< The englobing tetrahedron of each particle
    std::vector<bool> has_basis_function_;  ///< Marks active fem dofs vs inactive fem dofs

    std::vector<fem_interpolation_function_type>
        fem_interpolation_fields_; ///< Interpolations in interior tetrahedra
};

template <class KernelFunctionType>
inline fem_sph_model_t<KernelFunctionType>::fem_sph_model_t(
    geometry::tetrahedral_domain_t const& domain,
    geometry::grid_t const& grid,
    scalar_type support)
    : base_type(domain),
      sph_model_(),
      grid_(grid),
      Vjs_(),
      Wjs_(),
      Fjs_(),
      mixed_interpolation_fields_(),
      mixed_deformation_gradient_functions_(),
      particles_in_tet_(),
      particle_tets_(),
      has_basis_function_(),
      fem_interpolation_fields_()
{
    // Initialize dofs to material space positions
    has_basis_function_.resize(this->dof_count(), true);
    auto const& topology    = this->domain().topology();
    auto const vertex_count = topology.vertex_count();
    for (auto i = 0u; i < this->dof_count(); ++i)
    {
        this->dof(i) = this->point(i);
        // All points that coincide with tetrahedral mesh vertices
        // are stored first in the dof/point vector of tetrahedral_fem_model_t
        if (i < vertex_count)
        {
            bool const is_boundary_vertex = topology.is_boundary_vertex(i);
            if (is_boundary_vertex)
            {
                has_basis_function_[i] = false;
            }
        }
    }
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

    // Sample meshless nodes
    grid_.walk([this](Eigen::Vector3d const& X) {
        bool const is_grid_cell_center_in_domain = this->domain().contains(X);
        if (!is_grid_cell_center_in_domain)
            return;

        // Remove to test with particles everywhere in the domain
        index_type const ti = this->domain().in_tetrahedron(X);
        bool const is_in_boundary_tetrahedron =
            this->domain().topology().is_boundary_tetrahedron(ti);
        if (!is_in_boundary_tetrahedron)
            return;

        sph_model_.add_point(X);
        sph_model_.add_dof(X);
    });

    Eigen::Vector3d const& dX = grid_.delta();
    scalar_type const length  = dX.norm();
    scalar_type const h       = support * length;
    sph_model_.set_support_radius(h);
    sph_model_.initialize_in_support_query(
        length * 1e-2); // Add a 1% tolerance to the spatial search data structure based on the
                        // grid cells' dimensions

    std::vector<scalar_type> Vjs{};
    Vjs.reserve(sph_model_.point_count());
    std::vector<kernel_function_type> Wjs{};
    Wjs.reserve(sph_model_.point_count());
    std::vector<std::vector<index_type>> Njks{};
    Njks.reserve(sph_model_.point_count());

    particles_in_tet_.resize(this->domain().topology().tetrahedron_count());
    particle_tets_.reserve(sph_model_.point_count());

    // Create node shape functions and nodal neighbourhoods
    for (auto j = 0u; j < sph_model_.point_count(); ++j)
    {
        Eigen::Vector3d const& Xj = sph_model_.point(j);
        kernel_function_type Wj(autodiff::Vector3dual(Xj), h);
        Wjs.push_back(Wj);

        std::vector<index_type> nodes = sph_model_.in_support_of_nodes(Xj);
        Njks.push_back(nodes);

        index_type const ti = this->domain().in_tetrahedron(Xj);
        assert(ti != std::numeric_limits<index_type>::max());
        particles_in_tet_[ti].push_back(j);
        particle_tets_.push_back(ti);
    }
    for (auto j = 0u; j < sph_model_.point_count(); ++j)
    {
        Eigen::Vector3d const& Xj          = sph_model_.point(j);
        std::vector<index_type> const& Nks = Njks[j];
        scalar_type const density =
            std::accumulate(Nks.begin(), Nks.end(), 0., [&](scalar_type sum, index_type const Nk) {
                kernel_function_type const& Wk = Wjs[Nk];
                auto const dualW               = Wk(autodiff::Vector3dual(Xj));
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
    this->Fjs_.resize(sph_model_.point_count(), Eigen::Matrix3d::Identity());

    this->mixed_interpolation_fields_.reserve(sph_model_.point_count());
    this->mixed_deformation_gradient_functions_.reserve(sph_model_.point_count());

    for (auto j = 0u; j < sph_model_.point_count(); ++j)
    {
        Eigen::Vector3d const& Xj = sph_model_.point(j);
        index_type const ti       = particle_tets_[j];
        auto const& cell          = this->cell(ti);

        mixed_interpolation_function_type mixed_interpolation(
            Xj,
            ti,
            &(this->cells()),
            &(this->points()),
            &(this->dofs()),
            &(this->has_basis_function_),
            Njks[j],
            &(sph_model_.points()),
            &(sph_model_.dofs()),
            &(this->Vjs_),
            &(this->Wjs_),
            &(this->Fjs_));
        this->mixed_interpolation_fields_.push_back(mixed_interpolation);

        mixed_deformation_gradient_function_type mixed_deformation_gradient(
            j,
            this->mixed_interpolation_fields_.back());
        this->mixed_deformation_gradient_functions_.push_back(mixed_deformation_gradient);
    }
}

template <class KernelFunctionType>
inline fem_sph_model_t<KernelFunctionType>::fem_sph_model_t(self_type const& other)
    : base_type(other),
      sph_model_(other.sph_model_),
      grid_(other.grid_),
      Vjs_(other.Vjs_),
      Wjs_(other.Wjs_),
      Fjs_(other.Fjs_),
      mixed_interpolation_fields_(),
      mixed_deformation_gradient_functions_(),
      particles_in_tet_(other.particles_in_tet_),
      particle_tets_(other.particle_tets_),
      has_basis_function_(other.has_basis_function_),
      fem_interpolation_fields_()
{
    mixed_interpolation_fields_.reserve(other.sph_model_.point_count());
    mixed_deformation_gradient_functions_.reserve(other.sph_model_.point_count());
    for (auto j = 0u; j < other.sph_model_.point_count(); ++j)
    {
        Eigen::Vector3d const& Xj = other.point(j);
        mixed_interpolation_function_type const& other_interpolation =
            other.mixed_interpolation_fields_[j];
        index_type const ti = particle_tets_[j];
        mixed_interpolation_function_type interpolation(
            Xj,
            ti,
            &(this->cells()),
            &(this->points()),
            &(this->dofs()),
            &(this->has_basis_function_),
            other_interpolation.js,
            &(sph_model_.points()),
            &(sph_model_.dofs()),
            &(this->Vjs_),
            &(this->Wjs_),
            &(this->Fjs_));
        mixed_interpolation_fields_.push_back(interpolation);

        mixed_deformation_gradient_function_type deformation_gradient_function(
            j,
            this->mixed_interpolation_fields_[j]);
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
inline fem_sph_model_t<KernelFunctionType>&
fem_sph_model_t<KernelFunctionType>::operator=(self_type const& other)
{
    base_type::operator=(other);

    sph_model_          = other.sph_model_;
    grid_               = other.grid_;
    Vjs_                = other.Vjs_;
    Wjs_                = other.Wjs_;
    Fjs_                = other.Fjs_;
    particles_in_tet_   = other.particles_in_tet_;
    particle_tets_      = other.particle_tets_;
    has_basis_function_ = other.has_basis_function_;

    mixed_interpolation_fields_.reserve(other.sph_model_.point_count());
    mixed_deformation_gradient_functions_.reserve(other.sph_model_.point_count());
    for (auto j = 0u; j < other.sph_model_.point_count(); ++j)
    {
        Eigen::Vector3d const& Xj = other.sph_model_.point(j);
        mixed_interpolation_function_type const& other_interpolation =
            other.mixed_interpolation_fields_[j];
        index_type const ti = particle_tets_[j];
        auto const& cell    = this->cell(ti);
        mixed_interpolation_function_type interpolation(
            Xj,
            ti,
            &(this->cells()),
            &(this->points()),
            &(this->dofs()),
            &(this->has_basis_function_),
            other_interpolation.js,
            &(sph_model_.points()),
            &(sph_model_.dofs()),
            &(this->Vjs_),
            &(this->Wjs_),
            &(this->Fjs_));
        mixed_interpolation_fields_.push_back(interpolation);

        mixed_deformation_gradient_function_type deformation_gradient_function(
            j,
            this->mixed_interpolation_fields_[j]);
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
inline typename fem_sph_model_t<KernelFunctionType>::mixed_interpolation_function_type
fem_sph_model_t<KernelFunctionType>::mixed_interpolation_field_at(Eigen::Vector3d const& X)
{
    auto const& fem_domain             = this->domain();
    index_type const e                 = fem_domain.in_tetrahedron(X);
    std::vector<index_type> const& Njs = this->sph_model_.in_support_of_nodes(X);

    mixed_interpolation_function_type mixed_interpolation(
        X,
        e,
        &(this->cells()),
        &(this->points()),
        &(this->dofs()),
        &(this->has_basis_function_),
        Njs,
        &(sph_model_.points()),
        &(sph_model_.dofs()),
        &(this->Vjs_),
        &(this->Wjs_),
        &(this->Fjs_));

    return mixed_interpolation;
}

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_FEM_SPH_MODEL_H
