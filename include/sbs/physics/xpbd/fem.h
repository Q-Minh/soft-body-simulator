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
    math::green_strain_op_t strain_func_;
    math::stvk_strain_energy_density_op_t strain_energy_func_;
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
    Eigen::Matrix3d const E = strain_func_(F);
    scalar_type const Psi   = strain_energy_func_(E);

    scalar_type const C     = wg_ * Psi;
    Eigen::Matrix3d const P = strain_energy_func_.stress(F, E);

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

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_FEM_H