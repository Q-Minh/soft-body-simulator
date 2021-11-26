#ifndef SBS_PHYSICS_MECHANICS_SPH_MESHLESS_MODEL_H
#define SBS_PHYSICS_MECHANICS_SPH_MESHLESS_MODEL_H

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
class sph_meshless_model_t
    : public math::
          meshless_model_t<Eigen::Vector3d, Eigen::Vector3d, math::sph_basis_function_t<KernelType>>
{
  public:
    using kernel_function_type        = KernelType;
    using basis_function_type         = math::sph_basis_function_t<kernel_function_type>;
    using interpolation_function_type = math::sph_interpolation_t<kernel_function_type>;
    using deformation_gradient_function_type =
        math::sph_nodal_deformation_gradient_op_t<kernel_function_type>;
    using base_type = math::
        meshless_model_t<Eigen::Vector3d, Eigen::Vector3d, math::sph_basis_function_t<KernelType>>;
    using self_type = sph_meshless_model_t<kernel_function_type>;

    sph_meshless_model_t() = default;
    sph_meshless_model_t(
        geometry::tetrahedral_domain_t const& domain,
        geometry::grid_t const& grid,
        scalar_type support);

    sph_meshless_model_t(self_type const& other);
    sph_meshless_model_t& operator=(self_type const& other);

    // Accessors
    geometry::tetrahedral_domain_t const& domain() const { return domain_; }
    geometry::grid_t const& grid() const { return grid_; }

    scalar_type const& V(index_type i) const { return Vis_[i]; }
    kernel_function_type const& W(index_type i) const { return Wis_[i]; }
    Eigen::Matrix3d const& F(index_type i) const { return Fis_[i]; }
    Eigen::Matrix3d& F(index_type i) { return Fis_[i]; }

    interpolation_function_type const& interpolation_field_at(index_type i) const
    {
        return interpolation_fields_[i];
    }
    interpolation_function_type& interpolation_field_at(index_type i)
    {
        return interpolation_fields_[i];
    }
    interpolation_function_type interpolation_field_at(Eigen::Vector3d const& X) const;

    std::vector<index_type> const& neighbours(index_type i) const;

    deformation_gradient_function_type const& deformation_gradient_function(index_type i) const
    {
        return deformation_gradient_functions_[i];
    }

  private:
    geometry::tetrahedral_domain_t domain_; ///< The integration domain
    geometry::grid_t grid_;                 ///< The particle sampling grid
    std::vector<scalar_type> Vis_;          ///< Shepard coefficients (nodal volume)
    std::vector<kernel_function_type> Wis_; ///< Nodal kernel functions
    std::vector<Eigen::Matrix3d> Fis_;      ///< Deformation gradients at each point
    std::vector<interpolation_function_type>
        interpolation_fields_; ///< Interpolation fields at points
    std::vector<deformation_gradient_function_type>
        deformation_gradient_functions_; ///< Nodal deformation gradient callable functors
};

template <class KernelType>
inline sph_meshless_model_t<KernelType>::sph_meshless_model_t(
    geometry::tetrahedral_domain_t const& domain,
    geometry::grid_t const& grid,
    scalar_type support)
    : domain_(domain), grid_(grid), interpolation_fields_()
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

    std::vector<scalar_type> Vis{};
    Vis.reserve(this->point_count());
    std::vector<kernel_function_type> Wis{};
    Wis.reserve(this->point_count());
    std::vector<std::vector<index_type>> Nijs{};
    Nijs.reserve(this->point_count());

    // Create node shape functions and nodal neighbourhoods
    for (auto i = 0u; i < this->point_count(); ++i)
    {
        Eigen::Vector3d const& Xi = this->point(i);
        kernel_function_type Wi(autodiff::Vector3dual(Xi), h);
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
                auto const dualW               = Wj(autodiff::Vector3dual(Xi));
                scalar_type const W            = static_cast<scalar_type>(dualW);
                return sum + W;
            });
        scalar_type const Vi = 1. / density;
        Vis.push_back(Vi);
        kernel_function_type const& Wi = Wis[i];
        basis_function_type const sph_basis_function(Wi, Vi);
        this->add_basis_function(sph_basis_function);
    }

    this->Vis_ = Vis;
    this->Wis_ = Wis;
    this->Fis_.resize(this->point_count(), Eigen::Matrix3d::Identity());

    this->interpolation_fields_.reserve(this->point_count());
    this->deformation_gradient_functions_.reserve(this->point_count());
    for (auto i = 0u; i < this->point_count(); ++i)
    {
        Eigen::Vector3d const& Xi = this->point(i);
        interpolation_function_type sph_interpolation(
            Xi,
            Nijs[i],
            &(this->points()),
            &(this->dofs()),
            &(this->Vis_),
            &(this->Wis_),
            &(this->Fis_));

        this->interpolation_fields_.push_back(sph_interpolation);
        deformation_gradient_function_type sph_nodal_deformation_gradient(
            i,
            Xi,
            this->interpolation_fields_[i]);
        this->deformation_gradient_functions_.push_back(sph_nodal_deformation_gradient);
    }
}

template <class KernelType>
inline sph_meshless_model_t<KernelType>::sph_meshless_model_t(self_type const& other)
    : base_type(other),
      domain_(other.domain_),
      grid_(other.grid_),
      Vis_(other.Vis_),
      Wis_(other.Wis_),
      Fis_(other.Fis_),
      interpolation_fields_(),
      deformation_gradient_functions_()
{
    interpolation_fields_.reserve(other.point_count());
    deformation_gradient_functions_.reserve(other.point_count());
    for (auto i = 0u; i < other.point_count(); ++i)
    {
        Eigen::Vector3d const& Xi                              = other.point(i);
        interpolation_function_type const& other_interpolation = other.interpolation_fields_[i];
        interpolation_function_type interpolation(
            Xi,
            other_interpolation.js,
            &(this->points()),
            &(this->dofs()),
            &(this->Vis_),
            &(this->Wis_),
            &(this->Fis_));
        interpolation_fields_.push_back(interpolation);

        deformation_gradient_function_type deformation_gradient_function(
            i,
            Xi,
            this->interpolation_fields_[i]);
        this->deformation_gradient_functions_.push_back(deformation_gradient_function);
    }
}

template <class KernelType>
inline sph_meshless_model_t<KernelType>&
sph_meshless_model_t<KernelType>::operator=(self_type const& other)
{
    base_type::operator=(other);

    domain_ = other.domain_;
    grid_   = other.grid_;
    Vis_    = other.Vis_;
    Wis_    = other.Wis_;
    Fis_    = other.Fis_;

    interpolation_fields_.reserve(other.point_count());
    deformation_gradient_functions_.reserve(other.point_count());
    for (auto i = 0u; i < other.point_count(); ++i)
    {
        Eigen::Vector3d const& Xi                              = other.point(i);
        interpolation_function_type const& other_interpolation = other.interpolation_fields_[i];
        interpolation_function_type interpolation(
            Xi,
            other_interpolation.js,
            &(this->points()),
            &(this->dofs()),
            &(this->Vis_),
            &(this->Wis_),
            &(this->Fis_));
        interpolation_fields_.push_back(interpolation);

        deformation_gradient_function_type deformation_gradient_function(
            i,
            Xi,
            this->interpolation_fields_[i]);
        this->deformation_gradient_functions_.push_back(deformation_gradient_function);
    }
    return *this;
}

template <class KernelType>
inline typename sph_meshless_model_t<KernelType>::interpolation_function_type
sph_meshless_model_t<KernelType>::interpolation_field_at(Eigen::Vector3d const& X) const
{
    Eigen::Vector3d const Xi            = X;
    std::vector<index_type> const nodes = this->in_support_of_nodes(X);
    interpolation_function_type sph_interpolation(
        Xi,
        nodes,
        &(this->points()),
        &(this->dofs()),
        &(this->Vis_),
        &(this->Wis_),
        &(this->Fis_));
    return sph_interpolation;
}

template <class KernelType>
inline std::vector<index_type> const&
sph_meshless_model_t<KernelType>::neighbours(index_type i) const
{
    interpolation_function_type const& sph_interpolation = interpolation_field_at(i);
    return sph_interpolation.js;
}

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_SPH_MESHLESS_MODEL_H