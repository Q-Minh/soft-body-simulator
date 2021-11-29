#ifndef SBS_MATH_FEM_SPH_H
#define SBS_MATH_FEM_SPH_H

#include "sbs/math/kernels.h"

#include <Eigen/Core>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <vector>

namespace sbs {
namespace math {

template <class FemCellType, class KernelFunctionType = poly6_kernel_t>
struct fem_sph_interpolation_t
{
    using cell_type            = FemCellType;
    using kernel_function_type = KernelFunctionType;
    using self_type            = fem_sph_interpolation_t<FemCellType, KernelFunctionType>;

    fem_sph_interpolation_t() = default;
    fem_sph_interpolation_t(
        Eigen::Vector3d const& Xk,
        // FEM parameters
        index_type e,
        std::vector<cell_type> const* cells,
        std::vector<Eigen::Vector3d> const* Xis,
        std::vector<Eigen::Vector3d>* xis,
        std::vector<bool> const* has_basis_function,
        // SPH parameters
        std::vector<index_type> const& js,
        std::vector<Eigen::Vector3d> const* Xjs,
        std::vector<Eigen::Vector3d>* xjs,
        std::vector<scalar_type> const* Vjs,
        std::vector<kernel_function_type> const* Wjs,
        std::vector<Eigen::Matrix3d>* Fjs)
        : sk(0.),
          Xk(Xk),
          e(e),
          cells(cells),
          Xis(Xis),
          xis(xis),
          has_basis_function(has_basis_function),
          js(js),
          Xjs(Xjs),
          xjs(xjs),
          Vjs(Vjs),
          Wjs(Wjs),
          Fjs(Fjs)
    {
        sk = compute_shepard_coefficient(Xk);
        compute_is();
    }

    fem_sph_interpolation_t(self_type const& other) = default;
    fem_sph_interpolation_t& operator=(self_type const& other) = default;

    Eigen::Vector3d operator()(Eigen::Vector3d const& X) const
    {
        if (Xk.isApprox(X, sbs::eps()))
        {
            return eval();
        }
        else
        {
            return eval(X);
        }
    }

    Eigen::Vector3d eval() const { return eval(Xk, sk); }

    Eigen::Vector3d eval(Eigen::Vector3d const& X) const
    {
        scalar_type const s = compute_shepard_coefficient(X);
        return eval(X, s);
    }

    Eigen::Vector3d eval(Eigen::Vector3d const& X, scalar_type s) const
    {
        auto const& cell = (*cells)[e];

        Eigen::Vector3d xk{0., 0., 0.};

        // sum_i x_i phi_i(X)
        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            auto const i       = cell.node(r);
            bool const has_phi = (*has_basis_function)[i];
            if (has_phi)
            {
                auto const& phi           = cell.phi(r);
                Eigen::Vector3d const& xi = (*xis)[i];
                xk += xi * phi(X);
            }
        }
        // sum_j (F_j (X_k - X_j) + x_j) phi_j(X)
        for (index_type const j : js)
        {
            scalar_type const& Vj          = (*Vjs)[j];
            kernel_function_type const& Wj = (*Wjs)[j];
            scalar_type const Wkj          = static_cast<scalar_type>(Wj(X));
            Eigen::Matrix3d const& Fj      = (*Fjs)[j];
            Eigen::Vector3d const& Xj      = (*Xjs)[j];
            Eigen::Vector3d const& xj      = (*xjs)[j];
            Eigen::Vector3d const Xkj      = X - Xj;
            xk += Vj * (Fj * Xkj + xj) * Wkj;
        }
        // xk = sk * (sum_i x_i phi_i(X) + sum_j (F_j (X_k - X_j) + x_j) phi_j(X))
        xk = s * xk;
        return xk;
    }

    /**
     * @brief
     * Normally, the gradient/jacobian of a 3d vector valued function w.r.t. to a 3d parameter
     * should yield a 3x3 matrix. However, in this particular case, the jacobian becomes the
     * identity matrix scaled by (sk * phi_i) or (sk * Vj * Wkj), since dxk/dxk = I, and we have
     * that xk = (sk * sum_i x_i phi_i) + (sk * sum_j Vj (Fj*Xkj + xj) Wkj). As such, we return only
     * the scalars sk * Vj * Wkj for all neighbour nodes.
     * @return
     */
    std::pair<std::vector<scalar_type>, std::vector<scalar_type>> dxdxk() const
    {
        std::vector<scalar_type> dxdxis{};
        dxdxis.reserve(is.size());
        std::vector<scalar_type> dxdxjs{};
        dxdxjs.reserve(js.size());

        auto const& cell = (*cells)[e];
        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            auto const i       = cell.node(r);
            bool const has_phi = (*has_basis_function)[i];
            if (has_phi)
            {
                auto const& phi  = cell.phi(r);
                auto const dxdxi = sk * phi(Xk);
                dxdxis.push_back(dxdxi);
            }
        }
        for (index_type const j : js)
        {
            scalar_type const& Vj          = (*Vjs)[j];
            kernel_function_type const& Wj = (*Wjs)[j];
            scalar_type const Wkj          = static_cast<scalar_type>(Wj(Xk));
            scalar_type const dxdxj        = sk * Vj * Wkj;
            dxdxjs.push_back(dxdxj);
        }
        return std::make_pair(dxdxis, dxdxjs);
    }

    scalar_type compute_shepard_coefficient(Eigen::Vector3d const& X) const
    {
        scalar_type s = 0.;

        // sum_i phi_i(X)
        auto const& cell = (*cells)[e];
        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            auto const i       = cell.node(r);
            bool const has_phi = (*has_basis_function)[i];
            if (has_phi)
            {
                auto const& phi = cell.phi(r);
                s += phi(X);
            }
        }
        // sum_j phi_j(X)
        for (index_type const j : js)
        {
            scalar_type const& Vj          = (*Vjs)[j];
            kernel_function_type const& Wj = (*Wjs)[j];
            scalar_type const Wkj          = static_cast<scalar_type>(Wj(X));
            s += Vj * Wkj;
        }
        s = (1. / s);
        return s;
    }

    void compute_is()
    {
        is.clear();

        auto const& cell = (*cells)[e];
        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            auto const i       = cell.node(r);
            bool const has_phi = (*has_basis_function)[i];
            if (has_phi)
            {
                is.push_back(i);
            }
        }
    }

    scalar_type
        sk; ///< Shepard coefficient for maintaing 0-th order reproducibility of this interpolation
    Eigen::Vector3d Xk;                      ///< Point at which we evaluate this interpolation
    index_type e;                            ///< Index of cell in which this interpolation exists
    std::vector<cell_type> const* cells;     ///< The cells of the fem model
    std::vector<index_type> is;              ///< Indices of the fem basis functions that are active
    std::vector<Eigen::Vector3d> const* Xis; ///< Points of the mesh model
    std::vector<Eigen::Vector3d>* xis;       ///< Dofs of the mesh model
    std::vector<bool> const* has_basis_function;  ///< Flags for if a dof is active or not
    std::vector<index_type> js;                   ///< indices of meshless neighbours
    std::vector<Eigen::Vector3d> const* Xjs;      ///< Points of the meshless model
    std::vector<Eigen::Vector3d>* xjs;            ///< Dofs of the meshless model
    std::vector<scalar_type> const* Vjs;          ///< Volumes of the meshless model
    std::vector<kernel_function_type> const* Wjs; ///< Kernels of the meshless nodes
    std::vector<Eigen::Matrix3d>* Fjs;            ///< Deformation gradients at the meshless nodes
};

template <class FemCellType, class KernelFunctionType = poly6_kernel_t>
struct fem_sph_nodal_deformation_gradient_op_t
{
    using cell_type             = FemCellType;
    using kernel_function_type  = KernelFunctionType;
    using interpolation_op_type = fem_sph_interpolation_t<cell_type, kernel_function_type>;

    fem_sph_nodal_deformation_gradient_op_t(
        index_type k,
        interpolation_op_type const& interpolation_op)
        : idx_of_k(), k(k), interpolation_op(interpolation_op), gradWkjs(), gradphis(), Lk()
    {
        auto begin = interpolation_op.js.begin();
        auto end   = interpolation_op.js.end();
        auto it    = std::find(begin, end, k);
        idx_of_k   = static_cast<index_type>(std::distance(begin, it));

        // Precompute correction matrix L and basis function gradients
        store_L_and_grads();
    }

    void store_L_and_grads()
    {
        gradWkjs.clear();
        gradphis.clear();

        Lfem.setZero();
        Lsph.setZero();
        Lsphinv.setZero();
        Lk.setZero();

        auto const& Xis                = *interpolation_op.Xis;
        auto const& Xjs                = *interpolation_op.Xjs;
        auto const& cells              = *interpolation_op.cells;
        auto const& cell               = cells[interpolation_op.e];
        auto const& has_basis_function = *interpolation_op.has_basis_function;

        Eigen::Vector3d const& Xk = Xjs[k];

        // sum_i \nabla phi_i(X)
        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            auto const i       = cell.node(r);
            bool const has_phi = has_basis_function[i];
            if (has_phi)
            {
                auto const& phi               = cell.phi(r);
                Eigen::Vector3d const gradphi = phi.grad(Xk);
                gradphis.push_back(gradphi);

                Eigen::Vector3d const& Xi = Xis[i];
                Eigen::Vector3d const Xik = Xi - Xk;
                Lfem += gradphi * Xik.transpose();
            }
        }
        assert(interpolation_op.is.size() == gradphis.size());
        Eigen::Matrix3d const I = Eigen::Matrix3d::Identity();
        Lfem                    = I - Lfem;

        std::vector<scalar_type> const& Vjs          = *interpolation_op.Vjs;
        std::vector<kernel_function_type> const& Wjs = *interpolation_op.Wjs;

        // sum_j \nabla phi_j(X)
        auto const& js = interpolation_op.js;
        for (index_type const j : js)
        {
            scalar_type const& Vj          = Vjs[j];
            kernel_function_type const& Wj = Wjs[j];
            scalar_type const Wkj          = static_cast<scalar_type>(Wj(Xk));
            Eigen::Vector3d const& Xj      = Xjs[j];

            Eigen::Vector3d gradWkj = Wj.grad(Xk);
            gradWkjs.push_back(gradWkj);

            Eigen::Vector3d const Xjk = (Xj - Xk);
            Lsph += Vj * gradWkj * Xjk.transpose();
        }
#ifdef _DEBUG
        scalar_type const det = Lsph.determinant();
        assert(std::abs(det) > sbs::eps());
#endif
        Lsphinv = Lsph.inverse();

        Lk = Lfem * Lsphinv;
    }

    Eigen::Matrix3d eval() const
    {
        std::vector<index_type> const& is       = interpolation_op.is;
        std::vector<index_type> const& js       = interpolation_op.js;
        std::vector<Eigen::Vector3d> const& xis = *interpolation_op.xis;
        std::vector<Eigen::Vector3d> const& xjs = *interpolation_op.xjs;
        std::vector<scalar_type> const& Vjs     = *interpolation_op.Vjs;

        Eigen::Vector3d const xk = xjs[k];
        Eigen::Matrix3d F{};
        F.setZero();

        for (auto a = 0u; a < is.size(); ++a)
        {
            index_type const i             = is[a];
            Eigen::Vector3d const& gradphi = gradphis[a];
            Eigen::Vector3d const& xi      = xis[i];
            Eigen::Vector3d const xik      = xi - xk;
            F += xik * gradphi.transpose();
        }
        for (auto a = 0u; a < js.size(); ++a)
        {
            index_type const j             = js[a];
            Eigen::Vector3d const& gradWkj = gradWkjs[a];
            scalar_type const& Vj          = Vjs[j];
            Eigen::Vector3d const& xj      = xjs[j];
            Eigen::Vector3d const xjk      = xj - xk;
            F += xjk * (Vj * Lk * gradWkj).transpose();
        }

        return F;
    }

    /**
     * @brief
     * Returns a pair of the gradients w.r.t. fem dofs and gradients w.r.t. to sph dofs.
     * The format of the return value is the same as in the sph_nodal_deformation_gradient_op_t
     * @return
     */
    std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> dFdx() const
    {
        std::vector<index_type> const& is   = interpolation_op.is;
        std::vector<index_type> const& js   = interpolation_op.js;
        std::vector<scalar_type> const& Vjs = *interpolation_op.Vjs;

        std::vector<Eigen::Vector3d> dFdxis{};
        dFdxis.resize(is.size(), Eigen::Vector3d{0., 0., 0.});
        std::vector<Eigen::Vector3d> dFdxjs{};
        dFdxjs.resize(js.size(), Eigen::Vector3d{0., 0., 0.});

        for (auto a = 0u; a < is.size(); ++a)
        {
            index_type const i             = is[a];
            Eigen::Vector3d const& gradphi = gradphis[a];
            dFdxis[a]                      = gradphi;
            dFdxjs[idx_of_k] -= gradphi;
        }
        for (auto a = 0u; a < js.size(); ++a)
        {
            index_type const j = js[a];

            if (j == k)
                continue;

            Eigen::Vector3d const& gradWkj = gradWkjs[a];
            scalar_type const& Vj          = Vjs[j];
            Eigen::Vector3d const grad     = Lk * Vj * gradWkj;

            dFdxjs[a] = grad;
            dFdxjs[idx_of_k] -= grad;
        }

        return std::make_pair(dFdxis, dFdxjs);
    }

    index_type idx_of_k; ///< Position of the index k in interpolation_op.js
    index_type k;        ///< Index of the sph node
    interpolation_op_type const& interpolation_op; ///< The mixed sph+fem interpolation
    std::vector<Eigen::Vector3d> gradWkjs;         ///< Cached sph kernel gradients
    std::vector<Eigen::Vector3d> gradphis;         ///< Cached fem basis function gradients
                                                   ///< in this interpolation
    Eigen::Matrix3d Lfem;    ///< 3x3 fem contribution to the correction matrix
    Eigen::Matrix3d Lsph;    ///< 3x3 sph contribution to the correction matrix
    Eigen::Matrix3d Lsphinv; ///< Inverse of Lsph
    Eigen::Matrix3d Lk;      ///< Lsph^-1 * Lfem
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_FEM_SPH_H
