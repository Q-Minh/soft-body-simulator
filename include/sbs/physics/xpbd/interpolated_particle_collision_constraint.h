#ifndef SBS_PHYSICS_XPBD_INTERPOLATED_PARTICLE_COLLISION_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_INTERPOLATED_PARTICLE_COLLISION_CONSTRAINT_H

#include "sbs/physics/xpbd/constraint.h"
#include "sbs/physics/xpbd/simulation.h"

#include <Eigen/Core>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace sbs {
namespace physics {
namespace xpbd {

template <class InterpolationFunctionType>
class interpolated_particle_collision_constraint_t : public constraint_t
{
  public:
    using interpolation_op_type = InterpolationFunctionType;

    interpolated_particle_collision_constraint_t(
        scalar_type alpha /*compliance*/,
        scalar_type beta /*damping*/,
        std::vector<index_type> const& js /*particle indices*/,
        std::vector<index_type> const& bis /*particle body indices*/,
        interpolation_op_type const&
            interpolation_op /*function for interpolating the surface particle's position*/,
        Eigen::Vector3d const& Xi /*evaluation point of the interpolation*/,
        Eigen::Vector3d const& p, /*contact point*/
        Eigen::Vector3d const& n /*surface normal of correction*/);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    interpolation_op_type interpolation_op_;
    autodiff::Vector3dual Xi_; ///< Evaluation point of the interpolation operator

    Eigen::Vector3d qs_; ///< Intersection point
    Eigen::Vector3d n_;  ///< Normal at intersection point
};

template <class InterpolationFunctionType>
inline interpolated_particle_collision_constraint_t<InterpolationFunctionType>::
    interpolated_particle_collision_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        std::vector<index_type> const& js,
        std::vector<index_type> const& bis,
        interpolation_op_type const& interpolation_op,
        Eigen::Vector3d const& Xi,
        Eigen::Vector3d const& p,
        Eigen::Vector3d const& n)
    : constraint_t(alpha, beta, js, bis),
      interpolation_op_(interpolation_op),
      Xi_(Xi),
      qs_(p),
      n_(n)
{
}

template <class InterpolationFunctionType>
inline void
interpolated_particle_collision_constraint_t<InterpolationFunctionType>::project_positions(
    simulation_t& simulation,
    scalar_type dt)
{
    assert(js_.size() == interpolation_op_.uis.size());
    assert(interpolation_op_.uis.size() == interpolation_op_.phis.size());

    auto const& particles = simulation.particles();

    // Initialize interpolation field with latest particle positions
    for (auto i = 0u; i < js_.size(); ++i)
    {
        index_type const j        = js_[i];
        index_type const b        = bis_[i];
        particle_t const& p       = particles[b][j];
        Eigen::Vector3d const& xi = p.xi();
        interpolation_op_.uis[i]  = xi;
    }

    autodiff::Vector3dual x;
    auto const evalC =
        [this, &x](autodiff::Vector3dual Xi, Eigen::Vector3d const& qs, Eigen::Vector3d const& n) {
            x                        = interpolation_op_(Xi);
            autodiff::dual const sdf = (x - qs).dot(n);
            return sdf;
        };

    scalar_type const C = static_cast<scalar_type>(evalC(Xi_, qs_, n_));
    // Conditional constraint, only project if there is penetration
    if (C >= 0.)
        return;

    std::vector<Eigen::Vector3d> gradC{};
    gradC.reserve(js_.size());
    for (auto i = 0u; i < js_.size(); ++i)
    {
        using autodiff::at;
        using autodiff::gradient;
        using autodiff::wrt;

        autodiff::Vector3dual const dCdxi =
            gradient(evalC, wrt(interpolation_op_.uis[i]), at(Xi_, qs_, n_));
        gradC.push_back(dCdxi.cast<scalar_type>());
    }

    this->project_positions_with_dampling(simulation, C, gradC, dt);
}

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_INTERPOLATED_PARTICLE_COLLISION_CONSTRAINT_H