#ifndef SBS_PHYSICS_XPBD_FEM_MIXED_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_FEM_MIXED_CONSTRAINT_H

#include "constraint.h"
#include "sbs/aliases.h"
#include "sbs/math/elasticity.h"
#include "simulation.h"

namespace sbs {
namespace physics {
namespace xpbd {

template <class FemMixedModelType>
class stvk_fem_mixed_nodal_integration_constraint_t : public constraint_t
{
  public:
    using fem_mixed_model_type = FemMixedModelType;
    using fem_model_type       = typename fem_mixed_model_type::fem_model_type;
    using meshless_model_type  = typename fem_mixed_model_type::meshless_model_type;
    using interpolation_field_type =
        typename fem_mixed_model_type::mixed_interpolation_function_type;

    stvk_fem_mixed_nodal_integration_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        index_type nj,
        index_type bj,
        scalar_type Vj,
        index_type e,
        fem_mixed_model_type& fem_mixed_model,
        scalar_type E,
        scalar_type nu,
        std::size_t fem_particle_offset,
        std::size_t meshless_particle_offset);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    index_type j_;                          ///< Meshless particle of mixed model
    index_type b_;                          ///< affected body in particles
    scalar_type V_;                         ///< Nodal quadrature weight (volume based)
    index_type e_;                          ///< Index of cell in which this meshless node resides
    fem_mixed_model_type& fem_mixed_model_; ///< The mixed model
    // math::green_strain_op_t strain_op_;     ///< Strain functor
    // math::stvk_strain_energy_density_op_t
    //     strain_energy_density_op_; ///< Strain energy density functor
    math::small_strain_tensor_op_t strain_op_; ///< Strain functor
    math::corotational_linear_elasticity_strain_energy_density_op_t
        strain_energy_density_op_; ///< Strain energy density functor

    std::size_t fem_particle_offset_;
    std::size_t meshless_particle_offset_;
};

template <class FemMixedModelType>
inline stvk_fem_mixed_nodal_integration_constraint_t<FemMixedModelType>::
    stvk_fem_mixed_nodal_integration_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        index_type nj,
        index_type bj,
        scalar_type Vj,
        index_type e,
        fem_mixed_model_type& fem_mixed_model,
        scalar_type E,
        scalar_type nu,
        std::size_t fem_particle_offset,
        std::size_t meshless_particle_offset)
    : constraint_t(alpha, beta),
      j_(nj),
      b_(bj),
      V_(Vj),
      e_(e),
      fem_mixed_model_(fem_mixed_model),
      strain_op_(),
      strain_energy_density_op_(E, nu),
      fem_particle_offset_(fem_particle_offset),
      meshless_particle_offset_(meshless_particle_offset)
{
    auto const& interpolation_field = fem_mixed_model_.mixed_interpolation_field_at(j_);
    std::vector<index_type> js{};
    js.reserve(interpolation_field.is.size() + interpolation_field.js.size());

    for (index_type const i : interpolation_field.is)
    {
        index_type const particle_idx = static_cast<index_type>(fem_particle_offset_) + i;
        js.push_back(particle_idx);
    }
    for (index_type const j : interpolation_field.js)
    {
        index_type const particle_idx = static_cast<index_type>(meshless_particle_offset_) + j;
        js.push_back(particle_idx);
    }

    std::vector<index_type> bis{};
    bis.resize(js.size(), b_);

    this->set_indices(js);
    this->set_bodies(bis);
}

template <class FemMixedModelType>
inline void stvk_fem_mixed_nodal_integration_constraint_t<FemMixedModelType>::project_positions(
    simulation_t& simulation,
    scalar_type dt)
{
    auto const& particles                   = simulation.particles()[b_];
    interpolation_field_type& interpolation = fem_mixed_model_.mixed_interpolation_field_at(j_);
    std::vector<Eigen::Vector3d>& xis       = *interpolation.xis;
    std::vector<Eigen::Vector3d>& xjs       = *interpolation.xjs;
    for (index_type const i : interpolation.is)
    {
        Eigen::Vector3d const& xi = particles[fem_particle_offset_ + i].xi();
        xis[i]                    = xi;
    }
    for (index_type const j : interpolation.js)
    {
        Eigen::Vector3d const& xj = particles[meshless_particle_offset_ + j].xi();
        xjs[j]                    = xj;
    }

    auto const& deformation_gradient_op = fem_mixed_model_.mixed_deformation_gradient_function(j_);
    Eigen::Matrix3d const F             = deformation_gradient_op.eval();
    // Eigen::Matrix3d const E             = strain_op_(F);
    auto const [R, S]       = strain_op_.get_RS(F);
    Eigen::Matrix3d const E = strain_op_(S);
    scalar_type const Psi   = strain_energy_density_op_(E);
    scalar_type const C     = V_ * Psi;

    // Eigen::Matrix3d const P = strain_energy_density_op_.stress(F, E);
    Eigen::Matrix3d const P     = strain_energy_density_op_.stress(R, F);
    auto const [dFdxis, dFdxjs] = deformation_gradient_op.dFdx();
    assert(interpolation.is.size() == dFdxis.size());
    assert(interpolation.js.size() == dFdxjs.size());

    std::vector<Eigen::Vector3d> gradC{};
    gradC.reserve(dFdxis.size() + dFdxjs.size());

    for (auto a = 0u; a < dFdxis.size(); ++a)
    {
        Eigen::Vector3d grad{0., 0., 0.};
        grad(0) = V_ * P.row(0).dot(dFdxis[a]);
        grad(1) = V_ * P.row(1).dot(dFdxis[a]);
        grad(2) = V_ * P.row(2).dot(dFdxis[a]);
        gradC.push_back(grad);
    }
    for (auto a = 0u; a < dFdxjs.size(); ++a)
    {
        Eigen::Vector3d grad{};
        grad(0) = V_ * P.row(0).dot(dFdxjs[a]);
        grad(1) = V_ * P.row(1).dot(dFdxjs[a]);
        grad(2) = V_ * P.row(2).dot(dFdxjs[a]);
        gradC.push_back(grad);
    }

    this->project_positions_with_dampling(simulation, C, gradC, dt);

    fem_mixed_model_.F(j_) = F;
}

template <class FemMixedVisualModelType>
class fem_mixed_collision_constraint_t : public constraint_t
{
  public:
    using fem_mixed_visual_model_type = FemMixedVisualModelType;

    fem_mixed_collision_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        index_type vi,
        index_type b,
        fem_mixed_visual_model_type& fem_mixed_visual_model,
        Eigen::Vector3d const& c,
        Eigen::Vector3d const& n,
        std::size_t fem_particle_offset,
        std::size_t meshless_particle_offset);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    index_type vi_;
    index_type b_;
    fem_mixed_visual_model_type& fem_mixed_visual_model_;
    Eigen::Vector3d c_;
    Eigen::Vector3d n_;
    std::size_t fem_particle_offset_;
    std::size_t meshless_particle_offset_;
};

template <class FemMixedVisualModelType>
inline fem_mixed_collision_constraint_t<FemMixedVisualModelType>::fem_mixed_collision_constraint_t(
    scalar_type alpha,
    scalar_type beta,
    index_type vi,
    index_type b,
    fem_mixed_visual_model_type& fem_mixed_visual_model,
    Eigen::Vector3d const& c,
    Eigen::Vector3d const& n,
    std::size_t fem_particle_offset,
    std::size_t meshless_particle_offset)
    : constraint_t(alpha, beta),
      vi_(vi),
      b_(b),
      fem_mixed_visual_model_(fem_mixed_visual_model),
      c_(c),
      n_(n),
      fem_particle_offset_(fem_particle_offset),
      meshless_particle_offset_(meshless_particle_offset)
{
    auto const& interpolate = fem_mixed_visual_model_.interpolation_operator(vi_);

    std::vector<index_type> js{};
    js.reserve(interpolate.js.size() + interpolate.is.size());

    for (index_type const i : interpolate.is)
    {
        index_type const particle_idx = static_cast<index_type>(fem_particle_offset_) + i;
        js.push_back(particle_idx);
    }
    for (index_type const j : interpolate.js)
    {
        index_type const particle_idx = static_cast<index_type>(meshless_particle_offset_) + j;
        js.push_back(particle_idx);
    }

    std::vector<index_type> bis{};
    bis.resize(js.size(), b_);
    this->set_indices(js);
    this->set_bodies(bis);
}

template <class FemMixedVisualModelType>
inline void fem_mixed_collision_constraint_t<FemMixedVisualModelType>::project_positions(
    simulation_t& simulation,
    scalar_type dt)
{
    auto const& particles = simulation.particles()[b_];
    auto& interpolate     = fem_mixed_visual_model_.interpolation_operator(vi_);

    for (index_type const i : interpolate.is)
    {
        index_type const particle_idx = static_cast<index_type>(fem_particle_offset_) + i;
        Eigen::Vector3d const& xi     = particles[particle_idx].xi();
        (*interpolate.xis)[i]         = xi;
    }
    for (index_type const j : interpolate.js)
    {
        index_type const particle_idx = static_cast<index_type>(meshless_particle_offset_) + j;
        Eigen::Vector3d const& xj     = particles[particle_idx].xi();
        (*interpolate.xjs)[j]         = xj;
    }

    Eigen::Vector3d const p = interpolate.eval();

    scalar_type const C = (p - c_).dot(n_);

    if (C >= 0.)
    {
        return;
    }

    auto const [dpdxis, dpdxjs] = interpolate.dxdxk();

    assert(dpdxis.size() == interpolate.is.size());
    assert(dpdxjs.size() == interpolate.js.size());

    std::vector<Eigen::Vector3d> gradC{};
    gradC.reserve(dpdxis.size() + dpdxjs.size());
    for (std::size_t a = 0u; a < dpdxis.size(); ++a)
    {
        scalar_type const dpdxi    = dpdxis[a];
        Eigen::Vector3d const grad = dpdxi * n_;
        gradC.push_back(grad);
    }
    for (std::size_t a = 0u; a < dpdxjs.size(); ++a)
    {
        scalar_type const dpdxj    = dpdxjs[a];
        Eigen::Vector3d const grad = dpdxj * n_;
        gradC.push_back(grad);
    }

    this->project_positions_with_dampling(simulation, C, gradC, dt);
}

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_FEM_MIXED_CONSTRAINT_H
