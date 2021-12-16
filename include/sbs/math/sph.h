#ifndef SBS_MATH_SPH_H
#define SBS_MATH_SPH_H

#include "sbs/math/kernels.h"

#include <Eigen/Core>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <vector>

namespace sbs {
namespace math {

template <class KernelFunctionType = math::poly6_kernel_t>
struct sph_interpolation_t
{
    using kernel_function_type = KernelFunctionType;

    sph_interpolation_t(
        Eigen::Vector3d const& Xk,
        std::vector<index_type> const& js,
        std::vector<Eigen::Vector3d> const* Xjs,
        std::vector<Eigen::Vector3d> const* xjs,
        std::vector<scalar_type> const* Vjs,
        std::vector<kernel_function_type> const* Wjs)
        : sk(0.), Xk(Xk), js(js), Xjs(Xjs), xjs(xjs), Vjs(Vjs), Wjs(Wjs)
    {
        sk = compute_shepard_coefficient(Xk);
    }

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

    scalar_type compute_shepard_coefficient(Eigen::Vector3d const& X) const
    {
        scalar_type s = 0.;
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

    Eigen::Vector3d eval(Eigen::Vector3d const& X) const
    {
        scalar_type const s = compute_shepard_coefficient(X);
        return eval(X, s);
    }

    Eigen::Vector3d eval(Eigen::Vector3d const& X, scalar_type s) const
    {
        Eigen::Vector3d x{0., 0., 0.};

        for (index_type const j : js)
        {
            scalar_type const& Vj          = (*Vjs)[j];
            Eigen::Vector3d const& Xj      = (*Xjs)[j];
            Eigen::Vector3d const& xj      = (*xjs)[j];
            kernel_function_type const& Wj = (*Wjs)[j];
            scalar_type const Wkj          = static_cast<scalar_type>(Wj(X));
            Eigen::Vector3d const Xkj      = X - Xj;
            x += Vj * xj * Wkj;
        }
        x = s * x;

        return x;
    }

    /**
     * @brief
     * Normally, the gradient/jacobian of a 3d vector valued function w.r.t. to a 3d parameter
     * should yield a 3x3 matrix. However, in this particular case, the jacobian becomes the
     * identity matrix scaled by sk * Vj * Wkj, since dxk/dxk = I, and we have that
     * xk = sk * sum_j Vj (Fj*Xkj + xj) Wkj.
     * As such, we return only the scalars sk * Vj * Wkj for all neighbour nodes.
     * @return
     */
    std::vector<scalar_type> dxdxk() const
    {
        std::vector<scalar_type> dxdxjs{};
        dxdxjs.reserve(js.size());
        for (index_type const j : js)
        {
            scalar_type const& Vj          = (*Vjs)[j];
            kernel_function_type const& Wj = (*Wjs)[j];
            scalar_type const Wkj          = static_cast<scalar_type>(Wj(Xk));
            scalar_type const dxdxj        = sk * Vj * Wkj;
            dxdxjs.push_back(dxdxj);
        }
        return dxdxjs;
    }

    scalar_type sk;
    Eigen::Vector3d Xk;
    std::vector<index_type> js;
    std::vector<scalar_type> const* Vjs;
    std::vector<Eigen::Vector3d> const* Xjs;
    std::vector<Eigen::Vector3d> const* xjs;
    std::vector<kernel_function_type> const* Wjs;
};

template <class KernelFunctionType = math::poly6_kernel_t>
struct sph_nodal_deformation_gradient_op_t
{
    using kernel_function_type      = KernelFunctionType;
    using sph_interpolation_op_type = sph_interpolation_t<kernel_function_type>;

    sph_nodal_deformation_gradient_op_t(
        index_type i,
        sph_interpolation_op_type const& sph_interpolation_op)
        : idx_of_i(), i(i), sph_interpolation_op(sph_interpolation_op), gradWkjs(), Li()
    {
        auto begin = sph_interpolation_op.js.begin();
        auto end   = sph_interpolation_op.js.end();
        auto it    = std::find(begin, end, i);
        idx_of_i   = static_cast<index_type>(std::distance(begin, it));

        // Compute correction matrix L
        store_L_and_gradWijs();
    }

    void store_L_and_gradWijs()
    {
        auto const& js = sph_interpolation_op.js;
        // Cache basis function gradients and correction matrix L
        gradWkjs.clear();
        gradWkjs.reserve(js.size());
        Li.setZero();

        for (index_type const j : js)
        {
            using autodiff::at;
            using autodiff::gradient;
            using autodiff::wrt;

            kernel_function_type const& Wj = (*sph_interpolation_op.Wjs)[j];
            Eigen::Vector3d const& Xj      = (*sph_interpolation_op.Xjs)[j];
            Eigen::Vector3d const& Xi      = (*sph_interpolation_op.Xjs)[i];
            scalar_type const& Vj          = (*sph_interpolation_op.Vjs)[j];
            autodiff::Vector3dual dualXi(Xi);

            Eigen::Vector3d gradWij{0., 0., 0.};
            // Sph weight functions have nan gradients when Xi == Xj
            gradWij = Wj.grad(Xi);
            gradWkjs.push_back(gradWij);
            Eigen::Vector3d const Xji = (Xj - Xi);
            Li += Vj * gradWij * Xji.transpose();
        }
#ifdef _DEBUG
        scalar_type const det = Li.determinant();
#endif
        Li = Li.inverse().eval();
    }

    Eigen::Matrix3d eval() const
    {
        Eigen::Vector3d const xi = (*sph_interpolation_op.xjs)[i];
        Eigen::Matrix3d F{};
        F.setZero();
        std::vector<index_type> const& js = sph_interpolation_op.js;
        for (auto a = 0u; a < js.size(); ++a)
        {
            index_type const j             = js[a];
            Eigen::Vector3d const xj       = (*sph_interpolation_op.xjs)[j];
            Eigen::Vector3d const xji      = xj - xi;
            scalar_type const Vj           = (*sph_interpolation_op.Vjs)[j];
            Eigen::Vector3d const& gradWij = gradWkjs[a];
            F += xji * (Vj * Li * gradWij).transpose();
        }

        return F;
    }

    /**
     * @brief
     * dFdxk mathematically returns a 3x3x3 tensor. However, we can exploit knowledge of
     * the values in the tensor to return dFdxk more efficiently.
     * The derivative of F w.r.t. to the ith component of xk returns a 3x3 matrix of
     * all zeroes except for the ith row of the matrix. For any component i, the non-zero
     * row is always the same.
     * Essentially, (dF/dxk(0)).row(0) == (dF/dxk(1)).row(1) == (dF/dxk(2)).row(2) .
     * Hence, we just return the nonzero rows of the dF/dxk and it is up to the user to
     * interpret the return value correctly.
     * @return
     */
    std::vector<Eigen::Vector3d> dFdx() const
    {
        std::vector<Eigen::Vector3d> dFdxjs{};
        std::vector<index_type> const& js = sph_interpolation_op.js;
        dFdxjs.resize(js.size(), Eigen::Vector3d{0., 0., 0.});
        for (auto a = 0u; a < js.size(); ++a)
        {
            index_type const j = js[a];
            if (j == i)
                continue;

            auto const Vj                  = (*sph_interpolation_op.Vjs)[j];
            Eigen::Vector3d const& gradWij = gradWkjs[a];
            Eigen::Vector3d const dFdxj    = Vj * Li * gradWij;
            dFdxjs[a]                      = dFdxj;
            dFdxjs[idx_of_i] -= dFdxj;
        }
        return dFdxjs;
    }

    index_type idx_of_i;
    index_type i;
    sph_interpolation_op_type const&
        sph_interpolation_op; ///< The interpolation from that this gradient is computed from
    std::vector<Eigen::Vector3d> gradWkjs; ///< Cached basis function gradients
    Eigen::Matrix3d Li;                    ///< Cached correction matrix
};

template <class KernelFunctionType = math::poly6_kernel_t>
struct sph_efg_deformation_gradient_op_t
{
    using kernel_function_type      = KernelFunctionType;
    using sph_interpolation_op_type = sph_interpolation_t<kernel_function_type>;

    sph_efg_deformation_gradient_op_t(
        index_type k,
        std::vector<Eigen::Vector3d> const* Xks,
        std::vector<Eigen::Vector3d>* xks,
        sph_interpolation_op_type const& sph_interpolation_op)
        : sph_interpolation_op(sph_interpolation_op), gradWkjs(), Lk(), k(k), Xks(Xks), xks(xks)
    {
        // Compute correction matrix L
        store_L_and_gradWijs();
    }

    void store_L_and_gradWijs()
    {
        auto const& js = sph_interpolation_op.js;
        // Cache basis function gradients and correction matrix L
        gradWkjs.clear();
        gradWkjs.reserve(js.size());
        Lk.setZero();

        for (index_type const j : js)
        {
            using autodiff::at;
            using autodiff::gradient;
            using autodiff::wrt;

            kernel_function_type const& Wj = (*sph_interpolation_op.Wjs)[j];
            Eigen::Vector3d const& Xj      = (*sph_interpolation_op.Xjs)[j];
            Eigen::Vector3d const& Xk      = (*Xks)[k];
            scalar_type const& Vj          = (*sph_interpolation_op.Vjs)[j];

            Eigen::Vector3d gradWkj{0., 0., 0.};
            // Sph weight functions have nan gradients when Xk == Xj
            gradWkj = Wj.grad(Xk);
            gradWkjs.push_back(gradWkj);
            Eigen::Vector3d const Xjk = (Xj - Xk);
            Lk += Vj * gradWkj * Xjk.transpose();
        }
#ifdef _DEBUG
        scalar_type const det = Lk.determinant();
#endif
        Lk = Lk.inverse().eval();
    }

    Eigen::Matrix3d eval() const
    {
        Eigen::Vector3d const xk = (*xks)[k];
        Eigen::Matrix3d F{};
        F.setZero();
        std::vector<index_type> const& js = sph_interpolation_op.js;
        for (auto a = 0u; a < js.size(); ++a)
        {
            index_type const j             = js[a];
            Eigen::Vector3d const xj       = (*sph_interpolation_op.xjs)[j];
            Eigen::Vector3d const xjk      = xj - xk;
            scalar_type const Vj           = (*sph_interpolation_op.Vjs)[j];
            Eigen::Vector3d const& gradWkj = gradWkjs[a];
            F += xjk * (Vj * Lk * gradWkj).transpose();
        }

        return F;
    }

    /**
     * @brief
     * dFdxk mathematically returns a 3x3x3 tensor. However, we can exploit knowledge of
     * the values in the tensor to return dFdxk more efficiently.
     * The derivative of F w.r.t. to the ith component of xk returns a 3x3 matrix of
     * all zeroes except for the ith row of the matrix. For any component i, the non-zero
     * row is always the same.
     * Essentially, (dF/dxk(0)).row(0) == (dF/dxk(1)).row(1) == (dF/dxk(2)).row(2) .
     * Hence, we just return the nonzero rows of the dF/dxk and it is up to the user to
     * interpret the return value correctly.
     * @return
     */
    std::vector<Eigen::Vector3d> dFdx() const
    {
        std::vector<Eigen::Vector3d> dFdxjs{};
        std::vector<index_type> const& js = sph_interpolation_op.js;
        dFdxjs.resize(js.size(), Eigen::Vector3d{0., 0., 0.});
        for (auto a = 0u; a < js.size(); ++a)
        {
            index_type const j             = js[a];
            auto const Vj                  = (*sph_interpolation_op.Vjs)[j];
            Eigen::Vector3d const& gradWkj = gradWkjs[a];
            Eigen::Vector3d const dFdxj    = Vj * Lk * gradWkj;
            dFdxjs[a]                      = dFdxj;
        }
        return dFdxjs;
    }

    sph_interpolation_op_type const&
        sph_interpolation_op; ///< The interpolation from that this gradient is computed from
    std::vector<Eigen::Vector3d> gradWkjs; ///< Cached basis function gradients
    Eigen::Matrix3d Lk;                    ///< Cached correction matrix

    index_type k;                            ///< Index of integration point
    std::vector<Eigen::Vector3d> const* Xks; ///< Integration points
    std::vector<Eigen::Vector3d>* xks;       ///< Integration point coefficients
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_SPH_H