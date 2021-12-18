#ifndef SBS_PHYSICS_XPBD_SPH_H
#define SBS_PHYSICS_XPBD_SPH_H

#include "constraint.h"
#include "sbs/math/elasticity.h"
#include "simulation.h"

#include <cassert>

namespace sbs {
namespace physics {
namespace xpbd {

template <class SphMeshlessModelType>
class sph_integration_constraint_t : public constraint_t
{
  public:
    using sph_meshless_model_type = SphMeshlessModelType;

    sph_integration_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        index_type k,
        index_type bi,
        scalar_type Vi,
        sph_meshless_model_type& sph_model,
        scalar_type E,
        scalar_type nu);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    index_type k_;                       ///< Index of integration point
    index_type b_;                       ///< Index of sph body
    scalar_type V_;                      ///< Quadrature weight (nodal volume)
    sph_meshless_model_type& sph_model_; ///< The sph meshless mechanical model
    // math::green_strain_op_t strain_op_;  ///< The strain operator
    // math::stvk_strain_energy_density_op_t
    //     strain_energy_density_op_; ///< The strain energy density operator
    math::small_strain_tensor_op_t strain_op_; ///< The strain operator
    math::corotational_linear_elasticity_strain_energy_density_op_t
        strain_energy_density_op_; ///< The strain energy density operator
};

template <class SphMeshlessModelType>
inline sph_integration_constraint_t<SphMeshlessModelType>::sph_integration_constraint_t(
    scalar_type alpha,
    scalar_type beta,
    index_type k,
    index_type bi,
    scalar_type Vi,
    sph_meshless_model_type& sph_model,
    scalar_type E,
    scalar_type nu)
    : constraint_t(alpha, beta),
      k_(k),
      b_(bi),
      V_(Vi),
      sph_model_(sph_model),
      strain_op_(),
      strain_energy_density_op_(E, nu)
{
    std::vector<index_type> const js = sph_model_.neighbours(k_);
    std::vector<index_type> bis{};
    bis.resize(js.size(), b_);
    this->set_indices(js);
    this->set_bodies(bis);
}

template <class SphMeshlessModelType>
inline void sph_integration_constraint_t<SphMeshlessModelType>::project_positions(
    simulation_t& simulation,
    scalar_type dt)
{
    auto const& js        = this->js();
    auto const& particles = simulation.particles()[b_];
    for (index_type const j : js)
    {
        Eigen::Vector3d const& xj = particles[j].xi();
        sph_model_.dof(j)         = xj;
    }

    auto const& deformation_gradient_op = sph_model_.deformation_gradient_function(k_);
    Eigen::Matrix3d const F             = deformation_gradient_op.eval();
    // Eigen::Matrix3d const E             = strain_op_(F);
    auto const [R, S]       = strain_op_.get_RS(F);
    Eigen::Matrix3d const E = strain_op_(S);
    scalar_type const Psi   = strain_energy_density_op_(E);
    scalar_type const C     = V_ * Psi;

    // Eigen::Matrix3d const P               = strain_energy_density_op_.stress(F, E);
    Eigen::Matrix3d const P                   = strain_energy_density_op_.stress(R, F);
    std::vector<Eigen::Vector3d> const dFdxjs = deformation_gradient_op.dFdx();
    assert(js.size() == dFdxjs.size());
    std::vector<Eigen::Vector3d> gradC{};
    gradC.reserve(dFdxjs.size());

    for (auto a = 0u; a < dFdxjs.size(); ++a)
    {
        Eigen::Vector3d grad{};
        grad(0) = V_ * P.row(0).dot(dFdxjs[a]);
        grad(1) = V_ * P.row(1).dot(dFdxjs[a]);
        grad(2) = V_ * P.row(2).dot(dFdxjs[a]);
        gradC.push_back(grad);
    }

    this->project_positions_with_dampling(simulation, C, gradC, dt);
}

template <class SphVisualModelType>
class sph_collision_constraint_t : public constraint_t
{
  public:
    using sph_visual_model_type = SphVisualModelType;

    sph_collision_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        index_type vi,
        index_type b,
        sph_visual_model_type& sph_visual_model,
        Eigen::Vector3d const& c,
        Eigen::Vector3d const& n);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    index_type vi_;
    index_type b_;
    sph_visual_model_type& sph_visual_model_;
    Eigen::Vector3d c_;
    Eigen::Vector3d n_;
};

template <class SphVisualModelType>
inline sph_collision_constraint_t<SphVisualModelType>::sph_collision_constraint_t(
    scalar_type alpha,
    scalar_type beta,
    index_type vi,
    index_type b,
    sph_visual_model_type& sph_visual_model,
    Eigen::Vector3d const& c,
    Eigen::Vector3d const& n)
    : constraint_t(alpha, beta), vi_(vi), b_(b), sph_visual_model_(sph_visual_model), c_(c), n_(n)
{
    std::vector<index_type> const js = sph_visual_model_.meshless_neighbours_of_vertex(vi_);
    std::vector<index_type> bis{};
    bis.resize(js.size(), b_);
    this->set_indices(js);
    this->set_bodies(bis);
}

template <class SphVisualModelType>
inline void sph_collision_constraint_t<SphVisualModelType>::project_positions(
    simulation_t& simulation,
    scalar_type dt)
{
    auto const& js          = this->js();
    auto const& interpolate = sph_visual_model_.interpolation_operator(vi_);
    Eigen::Vector3d const p = interpolate.eval();

    scalar_type const C = (p - c_).dot(n_);

    if (C >= 0.)
    {
        return;
    }

    std::vector<scalar_type> dpdxk = interpolate.dxdxk();
    assert(dpdxk.size() == js.size());

    std::vector<Eigen::Vector3d> gradC{};
    gradC.reserve(js.size());
    for (std::size_t a = 0u; a < js.size(); ++a)
    {
        Eigen::Vector3d const grad = dpdxk[a] * n_;
        gradC.push_back(grad);
    }

    this->project_positions_with_dampling(simulation, C, gradC, dt);

    auto const& particles     = simulation.particles()[b_];
    auto mechanical_model_ptr = sph_visual_model_.mechanical_model();
    for (index_type const j : js)
    {
        Eigen::Vector3d const& xj    = particles[j].xi();
        mechanical_model_ptr->dof(j) = xj;
    }
}

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_SPH_H