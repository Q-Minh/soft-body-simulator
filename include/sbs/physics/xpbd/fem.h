#ifndef SBS_PHYSICS_XPBD_FEM_H
#define SBS_PHYSICS_XPBD_FEM_H

#include "constraint.h"
#include "sbs/aliases.h"
#include "sbs/math/elasticity.h"
#include "sbs/math/fem.h"
#include "simulation.h"

#include <Eigen/Core>
#include <vector>

namespace sbs {
namespace physics {
namespace xpbd {

template <class CellType>
class stvk_tetrahedral_quadrature_strain_constraint_t : public constraint_t
{
  public:
    using cell_type                          = CellType;
    using interpolation_function_type        = math::fem_interpolation_t<cell_type>;
    using deformation_gradient_function_type = math::fem_deformation_gradient_t<cell_type>;

    stvk_tetrahedral_quadrature_strain_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        std::size_t particle_offset,
        std::vector<index_type> const& is,
        std::vector<index_type> bis,
        interpolation_function_type const& interpolation_function,
        scalar_type wg,
        Eigen::Vector3d const& Xg,
        scalar_type E,
        scalar_type nu);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    std::size_t particle_offset_;
    interpolation_function_type interpolation_function_;
    scalar_type wg_;
    Eigen::Vector3d Xg_;
    // math::green_strain_op_t strain_func_;
    // math::stvk_strain_energy_density_op_t strain_energy_func_;
    math::small_strain_tensor_op_t strain_func_;
    math::corotational_linear_elasticity_strain_energy_density_op_t strain_energy_func_;
};

template <class CellType>
inline stvk_tetrahedral_quadrature_strain_constraint_t<CellType>::
    stvk_tetrahedral_quadrature_strain_constraint_t(
        scalar_type alpha,
        scalar_type beta,
        std::size_t particle_offset,
        std::vector<index_type> const& is,
        std::vector<index_type> bis,
        interpolation_function_type const& interpolation_function,
        scalar_type wg,
        Eigen::Vector3d const& Xg,
        scalar_type E,
        scalar_type nu)
    : constraint_t(alpha, beta),
      particle_offset_(particle_offset),
      interpolation_function_(interpolation_function),
      wg_(wg),
      Xg_(Xg),
      strain_func_(),
      strain_energy_func_(E, nu)
{
    std::vector<index_type> js{};
    for (index_type const i : is)
    {
        js.push_back(static_cast<index_type>(particle_offset) + i);
    }
    this->set_indices(js);
    this->set_bodies(bis);
}

template <class CellType>
inline void stvk_tetrahedral_quadrature_strain_constraint_t<CellType>::project_positions(
    simulation_t& simulation,
    scalar_type dt)
{
    auto const e      = interpolation_function_.e;
    auto const& cells = *interpolation_function_.cells;
    auto const& cell  = cells[e];
    assert(this->js().size() == cell.node_count());

    auto const& particles = simulation.particles();
    for (auto r = 0u; r < cell.node_count(); ++r)
    {
        index_type const i                = cell.node(r);
        index_type const b                = this->bs()[r];
        index_type const idx              = this->js()[r];
        particle_t const& p               = particles[b][idx];
        Eigen::Vector3d const& xi         = p.xi();
        (*interpolation_function_.xis)[i] = xi;
    }

    deformation_gradient_function_type deformation_gradient_func(interpolation_function_);
    Eigen::Matrix3d const F = deformation_gradient_func.eval(Xg_);
    auto const [R, S]       = strain_func_.get_RS(F);
    Eigen::Matrix3d const E = strain_func_(S);
    // Eigen::Matrix3d const E = strain_func_(F);
    scalar_type const Psi = strain_energy_func_(E);
    // scalar_type const Psi   = strain_energy_func_(E);

    scalar_type const C = wg_ * Psi;
    // Eigen::Matrix3d const P = strain_energy_func_.stress(F, E);
    Eigen::Matrix3d const P = strain_energy_func_.stress(R, F);

    std::vector<Eigen::Vector3d> const dFdxis = deformation_gradient_func.dFdx(Xg_);
    assert(dFdxis.size() == this->js().size());

    std::vector<Eigen::Vector3d> gradC{};
    gradC.reserve(this->js().size());

    for (auto a = 0u; a < this->js().size(); ++a)
    {
        index_type const i = this->js()[a];
        Eigen::Vector3d grad{};
        grad(0) = wg_ * P.row(0).dot(dFdxis[a]);
        grad(1) = wg_ * P.row(1).dot(dFdxis[a]);
        grad(2) = wg_ * P.row(2).dot(dFdxis[a]);
        gradC.push_back(grad);
    }

    this->project_positions_with_dampling(simulation, C, gradC, dt);
}

template <class InterpolationFunctionType>
class fem_collision_constraint_t : public constraint_t
{
  public:
    using interpolation_op_type = InterpolationFunctionType;

    fem_collision_constraint_t(
        scalar_type alpha /*compliance*/,
        scalar_type beta /*damping*/,
        std::vector<index_type> const& js /*particle indices*/,
        index_type b /*body idx*/,
        interpolation_op_type const&
            interpolation_op /*function for interpolating the surface particle's position*/,
        Eigen::Vector3d const& Xi /*evaluation point of the interpolation*/,
        Eigen::Vector3d const& p, /*contact point*/
        Eigen::Vector3d const& n /*surface normal of correction*/);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  private:
    interpolation_op_type interpolation_op_;
    index_type b_;
    Eigen::Vector3d Xi_; ///< Evaluation point of the interpolation operator

    Eigen::Vector3d qs_; ///< Intersection point
    Eigen::Vector3d n_;  ///< Normal at intersection point
};

template <class InterpolationFunctionType>
inline fem_collision_constraint_t<InterpolationFunctionType>::fem_collision_constraint_t(
    scalar_type alpha,
    scalar_type beta,
    std::vector<index_type> const& js,
    index_type b /*body idx*/,
    interpolation_op_type const& interpolation_op,
    Eigen::Vector3d const& Xi,
    Eigen::Vector3d const& p,
    Eigen::Vector3d const& n)
    : constraint_t(alpha, beta), interpolation_op_(interpolation_op), b_(b), Xi_(Xi), qs_(p), n_(n)
{
    this->set_indices(js);
    std::vector<index_type> bis{};
    bis.resize(js.size(), b_);
    this->set_bodies(bis);
}

template <class InterpolationFunctionType>
inline void fem_collision_constraint_t<InterpolationFunctionType>::project_positions(
    simulation_t& simulation,
    scalar_type dt)
{
    auto const& particles = simulation.particles();

    // Initialize interpolation field with latest particle positions
    auto const e     = interpolation_op_.e;
    auto const& cell = (*interpolation_op_.cells)[e];
    for (auto r = 0u; r < cell.node_count(); ++r)
    {
        auto const i                = cell.node(r);
        particle_t const& p         = particles[b_][i];
        Eigen::Vector3d const& xi   = p.xi();
        (*interpolation_op_.xis)[i] = xi;
    }

    Eigen::Vector3d const x = interpolation_op_.eval(Xi_);

    scalar_type const C = (x - qs_).dot(n_);

    if (C >= 0.)
        return;

    std::vector<scalar_type> const dxdxk = interpolation_op_.dxdxk(Xi_);
    assert(dxdxk.size() == cell.node_count());

    std::vector<Eigen::Vector3d> gradC{};
    gradC.reserve(cell.node_count());

    for (auto r = 0u; r < cell.node_count(); ++r)
    {
        Eigen::Vector3d const grad = dxdxk[r] * n_;
        gradC.push_back(grad);
    }

    this->project_positions_with_dampling(simulation, C, gradC, dt);
}

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_FEM_H