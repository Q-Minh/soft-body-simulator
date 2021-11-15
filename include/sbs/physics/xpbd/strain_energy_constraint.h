#ifndef SBS_PHYSICS_XPBD_STRAIN_ENERGY_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_STRAIN_ENERGY_CONSTRAINT_H

#include "constraint.h"
#include "sbs/aliases.h"
#include "sbs/math/elasticity.h"
#include "sbs/math/interpolation.h"
#include "sbs/physics/xpbd/simulation.h"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <vector>

namespace sbs {
namespace physics {
namespace xpbd {

template <class InterpolationFunctionType>
class strain_energy_quadrature_constraint_t : public constraint_t
{
  public:
    using interpolation_op_type = InterpolationFunctionType;

    strain_energy_quadrature_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        std::vector<index_type> const& js,
        std::vector<index_type> const& bis,
        interpolation_op_type const& interpolation_op,
        autodiff::dual const& wi,
        autodiff::Vector3dual const& Xi,
        scalar_type young_modulus,
        scalar_type poisson_ratio);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    interpolation_op_type interpolation_op_; ///< Interpolation field
    autodiff::dual wi_;                      ///< Quadrature weight
    autodiff::Vector3dual Xi_;               ///< Quadrature point
    scalar_type young_modulus_;
    scalar_type poisson_ratio_;
};

template <class InterpolationOpType>
inline strain_energy_quadrature_constraint_t<InterpolationOpType>::
    strain_energy_quadrature_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        std::vector<index_type> const& js,
        std::vector<index_type> const& bis,
        interpolation_op_type const& interpolation_op,
        autodiff::dual const& wi,
        autodiff::Vector3dual const& Xi,
        scalar_type young_modulus,
        scalar_type poisson_ratio)
    : constraint_t(alpha, beta, js, bis),
      interpolation_op_(interpolation_op),
      wi_(wi),
      Xi_(Xi),
      young_modulus_(young_modulus),
      poisson_ratio_(poisson_ratio)
{
}

template <class InterpolationOpType>
inline void strain_energy_quadrature_constraint_t<InterpolationOpType>::project_positions(
    simulation_t& simulation,
    scalar_type dt)
{
    assert(js_.size() == interpolation_op_.uis.size());
    assert(interpolation_op_.uis.size() == interpolation_op_.phis.size());

    // Initialize interpolation field with current particle positions
    auto const& particles = simulation.particles();
    for (auto i = 0u; i < js_.size(); ++i)
    {
        index_type const j        = js_[i];
        index_type const b        = bis_[i];
        particle_t const& p       = particles[b][j];
        Eigen::Vector3d const& xi = p.xi();
        interpolation_op_.uis[i]  = xi;
    }

    using deformation_gradient_op_type =
        sbs::math::deformation_gradient_op_t<interpolation_op_type>;
    using strain_op_type = sbs::math::strain_op_t<deformation_gradient_op_type>;

    deformation_gradient_op_type deformation_gradient_op(this->interpolation_op_);
    strain_op_type strain_op(deformation_gradient_op);
    sbs::math::strain_energy_density_op_t<strain_op_type> strain_energy_density_op(
        strain_op,
        young_modulus_,
        poisson_ratio_);

    autodiff::Matrix3dual F, E;

    auto const total_energy = [&F, &E, &strain_energy_density_op](
                                  autodiff::dual wi,
                                  autodiff::Vector3dual Xi) -> autodiff::dual {
        return wi * strain_energy_density_op(Xi, F, E);
    };

    autodiff::dual wiPsi = total_energy(wi_, Xi_);
    // Evaluate constraint
    scalar_type const C = static_cast<scalar_type>(wiPsi);

    // Constraint projection
    std::vector<Eigen::Vector3d> gradC{};
    gradC.reserve(js_.size());

    using autodiff::at;
    using autodiff::gradient;
    using autodiff::wrt;

    for (auto i = 0u; i < js_.size(); ++i)
    {
        autodiff::Vector3dual& xi         = this->interpolation_op_.uis[i];
        autodiff::Vector3dual const dCdxi = autodiff::gradient(total_energy, wrt(xi), at(wi_, Xi_));
        gradC.push_back(dCdxi.cast<scalar_type>());
    }

    this->project_positions_with_dampling(simulation, C, gradC, dt);
}

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_STRAIN_ENERGY_CONSTRAINT_H