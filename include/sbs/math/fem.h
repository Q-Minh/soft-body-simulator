#ifndef SBS_MATH_FEM_H
#define SBS_MATH_FEM_H

#include "sbs/aliases.h"

#include <Eigen/Core>
#include <vector>

namespace sbs {
namespace math {

template <class CellType>
struct fem_interpolation_t
{
    using cell_type           = CellType;
    using basis_function_type = typename cell_type::basis_function_type;

    fem_interpolation_t(cell_type const& cell, std::vector<Eigen::Vector3d> const* xis)
        : cell(cell), xis(xis)
    {
    }

    Eigen::Vector3d eval(Eigen::Vector3d const& X) const
    {
        Eigen::Vector3d x{0., 0., 0.};
        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            index_type const i             = cell.node(r);
            basis_function_type const& phi = cell.phi(r);
            Eigen::Vector3d const& xi      = (*xis)[i];
            x += xi * phi(X);
        }
        return x;
    }

    cell_type const& cell;
    std::vector<Eigen::Vector3d> const* xis;
};

template <class CellType>
struct fem_deformation_gradient_t
{
    using cell_type           = CellType;
    using interpolation_type  = fem_interpolation_t<cell_type>;
    using basis_function_type = typename cell_type::basis_function_type;

    fem_deformation_gradient_t(interpolation_type const& interpolation_op)
        : interpolation_op(interpolation_op)
    {
    }

    Eigen::Matrix3d eval(Eigen::Vector3d const& X) const
    {
        Eigen::Matrix3d F{};
        F.setZero();

        cell_type const& cell                   = interpolation_op.cell;
        std::vector<Eigen::Vector3d> const& xis = *interpolation_op.xis;

        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            index_type const i             = cell.node(r);
            basis_function_type const& phi = cell.phi(r);
            Eigen::Vector3d const& xi      = xis[i];
            Eigen::Vector3d const gradphi  = phi.grad(X);
            F += xi * gradphi.transpose();
        }

        return F;
    }

    interpolation_type const& interpolation_op;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_FEM_H