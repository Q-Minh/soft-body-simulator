#ifndef SBS_PHYSICS_MECHANICS_SPH_EFG_MODEL_H
#define SBS_PHYSICS_MECHANICS_SPH_EFG_MODEL_H

#include "sbs/geometry/grid.h"
#include "sbs/geometry/tetrahedral_domain.h"
#include "sbs/math/basis_functions.h"
#include "sbs/math/meshless_model.h"
#include "sbs/math/sph.h"

#include <Eigen/Core>
#include <vector>

namespace sbs {
namespace physics {
namespace mechanics {

template <class KernelType>
class sph_efg_model_t : public math::meshless_model_t<
                            Eigen::Vector3d,
                            Eigen::Vector3d,
                            math::differentiable::sph_basis_function_t<KernelType>>
{
  public:
    using kernel_function_type = KernelType;

    using basis_function_type = math::differentiable::sph_basis_function_t<kernel_function_type>;

    using interpolation_function_type = math::sph_efg_interpolation_op_t<kernel_function_type>;

    using deformation_gradient_function_type =
        math::sph_efg_deformation_gradient_op_t<kernel_function_type>;

    using base_type = math::meshless_model_t<
        Eigen::Vector3d,
        Eigen::Vector3d,
        math::differentiable::sph_basis_function_t<KernelType>>;

    using self_type = sph_efg_model_t<kernel_function_type>;

    sph_efg_model_t() = default;
    sph_efg_model_t(
        geometry::tetrahedral_domain_t const& domain,
        geometry::grid_t const& grid,
        scalar_type support);

    sph_efg_model_t(self_type const& other);
    sph_efg_model_t& operator=(self_type const& other);

    // Accessors
    geometry::tetrahedral_domain_t const& domain() const { return domain_; }
    geometry::grid_t const& grid() const { return grid_; }

    scalar_type const& V(index_type i) const { return Vis_[i]; }
    kernel_function_type const& W(index_type i) const { return Wis_[i]; }

    std::size_t integration_point_count() const { return Xks_.size(); }
    Eigen::Vector3d const& integration_point(index_type k) const { return Xks_[k]; }

    interpolation_function_type const& interpolation_field_at(index_type k) const
    {
        return interpolation_fields_[k];
    }
    interpolation_function_type& interpolation_field_at(index_type k)
    {
        return interpolation_fields_[k];
    }
    interpolation_function_type interpolation_field_at(Eigen::Vector3d const& X) const;

    std::vector<index_type> const& neighbours(index_type k) const;

    deformation_gradient_function_type const& deformation_gradient_function(index_type k) const
    {
        return deformation_gradient_functions_[k];
    }

    index_type tetrahedron_of_integration_point(index_type const k) const
    {
        return efg_point_tets_[k];
    }

    // modifiers
    void add_integration_point(Eigen::Vector3d const& X) { Xks_.push_back(X); }

    /**
     * @brief Update our world space positions at the integration points.
     */
    void update();

  private:
    geometry::tetrahedral_domain_t domain_; ///< The integration domain
    std::vector<scalar_type> Vis_;          ///< Shepard coefficients (nodal volume)
    std::vector<kernel_function_type> Wis_; ///< Nodal kernel functions

    std::vector<Eigen::Vector3d> Xks_; ///< Integration points

    std::vector<interpolation_function_type>
        interpolation_fields_; ///< Interpolation fields at efg integration points
    std::vector<deformation_gradient_function_type>
        deformation_gradient_functions_; ///< Deformation gradient callable functors at efg
                                         ///< integration points

    std::vector<index_type> efg_point_tets_; ///< The englobing tetrahedron of each efg point
};

template <class KernelType>
inline sph_efg_model_t<KernelType>::sph_efg_model_t(
    geometry::tetrahedral_domain_t const& domain,
    geometry::grid_t const& grid,
    scalar_type support)
    : domain_(domain),
      Vis_(),
      Wis_(),
      Xks_(),
      interpolation_fields_(),
      deformation_gradient_functions_(),
      efg_point_tets_()
{
    // Sample meshless nodes
    auto const& topology    = this->domain().topology();
    auto const vertex_count = topology.vertex_count();
    for (auto i = 0u; i < vertex_count; ++i)
    {
        index_type const vi       = static_cast<index_type>(i);
        Eigen::Vector3d const& Xi = this->domain().position(vi);
        this->add_point(Xi);
        this->add_dof(Xi);
    }
    // Sample efg integration points
    auto const tet_count = topology.tetrahedron_count();
    efg_point_tets_.reserve(tet_count);
    Xks_.reserve(tet_count);
    for (auto i = 0u; i < tet_count; ++i)
    {
        index_type const ti              = static_cast<index_type>(i);
        topology::tetrahedron_t const& t = topology.tetrahedron(ti);
        Eigen::Vector3d const& X1        = this->domain().position(t.v1());
        Eigen::Vector3d const& X2        = this->domain().position(t.v2());
        Eigen::Vector3d const& X3        = this->domain().position(t.v3());
        Eigen::Vector3d const& X4        = this->domain().position(t.v4());

        // Get barycenter of tetrahedron
        Eigen::Vector3d const Xk = 0.25 * (X1 + X2 + X3 + X4);
        this->add_integration_point(Xk);
        efg_point_tets_.push_back(i);
    }
    // Compute support radius of shape functions using largest edge length
    // in the tetrahedralization
    auto const edge_count       = topology.edge_count();
    scalar_type max_edge_length = 0.;
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

    scalar_type const h = max_edge_length;
    this->set_support_radius(h);
    this->initialize_in_support_query(h * 1e-2); // Add a 1% tolerance to the spatial search data
                                                 // structure based on the grid cells' dimensions

    std::vector<scalar_type> Vis{};
    Vis.reserve(this->point_count());
    std::vector<kernel_function_type> Wis{};
    Wis.reserve(this->point_count());
    std::vector<std::vector<index_type>> Nijs{};
    Nijs.reserve(this->point_count());

    // Create nodal shape functions and nodal neighbourhoods
    for (auto i = 0u; i < this->point_count(); ++i)
    {
        Eigen::Vector3d const& Xi = this->point(i);
        kernel_function_type Wi(Xi, h);
        Wis.push_back(Wi);

        std::vector<index_type> nodes = this->in_support_of_nodes(Xi);
        Nijs.push_back(nodes);
    }
    for (auto i = 0u; i < this->point_count(); ++i)
    {
        Eigen::Vector3d const& Xi          = this->point(i);
        std::vector<index_type> const& Njs = Nijs[i];
        scalar_type const density =
            std::accumulate(Njs.begin(), Njs.end(), 0., [&](scalar_type sum, index_type const Nj) {
                kernel_function_type const& Wj = Wis[Nj];
                scalar_type const W            = Wj(Xi);
                return sum + W;
            });
        scalar_type const Vi = 1. / density;
        Vis.push_back(Vi);
        kernel_function_type const& Wi = Wis[i];
        basis_function_type const sph_basis_function(Wi, Vi);
        this->add_basis_function(sph_basis_function);
    }

    // Create interpolation fields and deformation gradient fields at efg points
    std::vector<std::vector<index_type>> Nkjs{};
    Nkjs.reserve(this->integration_point_count());
    for (auto k = 0u; k < this->integration_point_count(); ++k)
    {
        Eigen::Vector3d const& Xk     = this->integration_point(k);
        std::vector<index_type> nodes = this->in_support_of_nodes(Xk);
        Nkjs.push_back(nodes);
        assert(nodes.size() >= 3u);
    }

    this->Vis_ = Vis;
    this->Wis_ = Wis;

    this->interpolation_fields_.reserve(this->integration_point_count());
    this->deformation_gradient_functions_.reserve(this->integration_point_count());
    for (auto k = 0u; k < this->integration_point_count(); ++k)
    {
        Eigen::Vector3d const& Xk = this->integration_point(k);
        interpolation_function_type sph_interpolation(
            Xk,
            Nkjs[k],
            &(this->points()),
            &(this->dofs()),
            &(this->Vis_),
            &(this->Wis_));

        this->interpolation_fields_.push_back(sph_interpolation);
        deformation_gradient_function_type sph_efg_deformation_gradient(
            k,
            &(this->Xks_),
            this->interpolation_fields_[k]);
        this->deformation_gradient_functions_.push_back(sph_efg_deformation_gradient);
    }
}

template <class KernelType>
inline sph_efg_model_t<KernelType>::sph_efg_model_t(self_type const& other)
    : base_type(other),
      domain_(other.domain_),
      Vis_(other.Vis_),
      Wis_(other.Wis_),
      Xks_(other.Xks_),
      interpolation_fields_(),
      deformation_gradient_functions_(),
      efg_point_tets_(other.efg_point_tets_)
{
    interpolation_fields_.reserve(other.integration_point_count());
    deformation_gradient_functions_.reserve(other.integration_point_count());
    for (auto k = 0u; k < other.integration_point_count(); ++k)
    {
        Eigen::Vector3d const& Xk                              = other.integration_point(k);
        interpolation_function_type const& other_interpolation = other.interpolation_fields_[k];
        interpolation_function_type interpolation(
            Xk,
            other_interpolation.js,
            &(this->points()),
            &(this->dofs()),
            &(this->Vis_),
            &(this->Wis_));
        interpolation_fields_.push_back(interpolation);

        deformation_gradient_function_type deformation_gradient_function(
            k,
            &(this->Xks_),
            this->interpolation_fields_[k]);
        this->deformation_gradient_functions_.push_back(deformation_gradient_function);
    }
}

template <class KernelType>
inline sph_efg_model_t<KernelType>& sph_efg_model_t<KernelType>::operator=(self_type const& other)
{
    base_type::operator=(other);

    domain_         = other.domain_;
    Vis_            = other.Vis_;
    Wis_            = other.Wis_;
    Xks_            = other.Xks_;
    efg_point_tets_ = other.efg_point_tets_;

    interpolation_fields_.reserve(other.integration_point_count());
    deformation_gradient_functions_.reserve(other.integration_point_count());
    for (auto k = 0u; k < other.integration_point_count(); ++k)
    {
        Eigen::Vector3d const& Xk                              = other.integration_point(k);
        interpolation_function_type const& other_interpolation = other.interpolation_fields_[k];
        interpolation_function_type interpolation(
            Xk,
            other_interpolation.js,
            &(this->points()),
            &(this->dofs()),
            &(this->Vis_),
            &(this->Wis_));
        interpolation_fields_.push_back(interpolation);

        deformation_gradient_function_type deformation_gradient_function(
            k,
            &(this->Xks_),
            this->interpolation_fields_[k]);
        this->deformation_gradient_functions_.push_back(deformation_gradient_function);
    }
    return *this;
}

template <class KernelType>
inline typename sph_efg_model_t<KernelType>::interpolation_function_type
sph_efg_model_t<KernelType>::interpolation_field_at(Eigen::Vector3d const& X) const
{
    std::vector<index_type> const nodes = this->in_support_of_nodes(X);
    interpolation_function_type sph_interpolation(
        X,
        nodes,
        &(this->points()),
        &(this->dofs()),
        &(this->Vis_),
        &(this->Wis_));
    return sph_interpolation;
}

template <class KernelType>
inline std::vector<index_type> const& sph_efg_model_t<KernelType>::neighbours(index_type k) const
{
    interpolation_function_type const& sph_interpolation = interpolation_field_at(k);
    return sph_interpolation.js;
}

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_SPH_EFG_MODEL_H
