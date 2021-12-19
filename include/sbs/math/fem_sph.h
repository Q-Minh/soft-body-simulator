#ifndef SBS_MATH_FEM_SPH_H
#define SBS_MATH_FEM_SPH_H

#include "sbs/math/kernels.h"

#include <Eigen/Core>
#include <array>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <vector>

namespace sbs {
namespace math {

template <class FemCellType, class KernelFunctionType = poly6_kernel_t>
struct fem_sph_efg_interpolation_t
{
    using cell_type            = FemCellType;
    using kernel_function_type = KernelFunctionType;
    using self_type            = fem_sph_efg_interpolation_t<FemCellType, KernelFunctionType>;

    fem_sph_efg_interpolation_t() = default;
    fem_sph_efg_interpolation_t(
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
        Eigen::Matrix3d const* Fk)
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
          Fk(Fk)
    {
        sk = compute_shepard_coefficient(Xk);
        compute_is();
    }

    fem_sph_efg_interpolation_t(self_type const& other) = default;
    fem_sph_efg_interpolation_t& operator=(self_type const& other) = default;

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
        Eigen::Vector3d const u = get_u(X);
        return s * u;
    }

    /**
     * @brief
     * Normally, the gradient/jacobian of a 3d vector valued function w.r.t. to a 3d parameter
     * should yield a 3x3 matrix. However, in this particular case, the jacobian becomes the
     * identity matrix scaled by (sk * phi_i) or (sk * Vj * Wkj), since dxk/dxk = I, and we have
     * that xk = (sk * sum_i x_i phi_i) + (sk * sum_j Vj xj Wkj). As such, we return only
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

    /**
     * @brief
     * Let x(X) be the 0th-order reproducible version of this interpolation
     * field. We reformulate x(X) as x(X) = u(X) / v(X), where
     * u(X) = \sum_i x_i \phi_i (X) + \sum_j x_j \phi_j (X),
     * v(X) = \sum_i \phi_i (X) + \sum_j \phi_j (X).
     * This method returns u(X).
     *
     * @param X Material space position.
     * @return u(X)
     */
    Eigen::Vector3d get_u(Eigen::Vector3d const& X) const
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
        // sum_j (x_j) phi_j(X)
        for (index_type const j : js)
        {
            scalar_type const& Vj          = (*Vjs)[j];
            kernel_function_type const& Wj = (*Wjs)[j];
            scalar_type const Wkj          = static_cast<scalar_type>(Wj(X));
            Eigen::Vector3d const& xj      = (*xjs)[j];
            Eigen::Vector3d const& Xj      = (*Xjs)[j];
            xk += (xj + (*Fk) * (X - Xj)) * Vj * Wkj;
        }
        return xk;
    }

    /**
     * @brief
     * Let x(X) be the 0th-order reproducible version of this interpolation
     * field. We reformulate x(X) as x(X) = u(X) / v(X), where
     * u(X) = \sum_i x_i \phi_i (X) + \sum_j x_j \phi_j (X),
     * v(X) = \sum_i \phi_i (X) + \sum_j \phi_j (X).
     * This method returns v(X).
     *
     * @param X Material space position X
     * @return v(X)
     */
    scalar_type get_v(Eigen::Vector3d const& X) const
    {
        scalar_type v = 0.;

        // sum_i phi_i(X)
        auto const& cell = (*cells)[e];
        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            auto const i       = cell.node(r);
            bool const has_phi = (*has_basis_function)[i];
            if (has_phi)
            {
                auto const& phi = cell.phi(r);
                v += phi(X);
            }
        }
        // sum_j phi_j(X)
        for (index_type const j : js)
        {
            scalar_type const& Vj          = (*Vjs)[j];
            kernel_function_type const& Wj = (*Wjs)[j];
            scalar_type const Wkj          = static_cast<scalar_type>(Wj(X));
            v += Vj * Wkj;
        }
        return v;
    }

    scalar_type compute_shepard_coefficient(Eigen::Vector3d const& X) const
    {
        scalar_type s = get_v(X);
        s             = (1. / s);
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

    scalar_type sk;     ///< Shepard coefficient for maintaing 0-th order reproducibility of this
                        ///< interpolation
    Eigen::Vector3d Xk; ///< Point at which we evaluate this interpolation
    index_type e;       ///< Index of cell in which this interpolation exists
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
    Eigen::Matrix3d const* Fk;                    ///< Deformation gradient in cell
};

template <class FemCellType, class KernelFunctionType = poly6_kernel_t>
struct fem_sph_nodal_interpolation_t
    : public fem_sph_efg_interpolation_t<FemCellType, KernelFunctionType>
{
    using cell_type            = FemCellType;
    using kernel_function_type = KernelFunctionType;
    using base_type            = fem_sph_efg_interpolation_t<FemCellType, KernelFunctionType>;
    using self_type            = fem_sph_nodal_interpolation_t<FemCellType, KernelFunctionType>;

    fem_sph_nodal_interpolation_t() = default;
    fem_sph_nodal_interpolation_t(
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
        : base_type(Xk, e, cells, Xis, xis, has_basis_function, js, Xjs, xjs, Vjs, Wjs, nullptr),
          Fjs(Fjs)
    {
    }

    fem_sph_nodal_interpolation_t(self_type const& other) = default;
    fem_sph_nodal_interpolation_t& operator=(self_type const& other) = default;

    /**
     * @brief
     * Solves the optimization problem argmin_X f(X), where
     * f(X) = 1/2 || x(X) - c ||^2 using Newton's method.
     * An initial guess Xk must be provided, as well as the
     * target point in world space c.
     * Essentially, this method gives the inverse mapping x^-1(X).
     * Given the point c in world space, we attempt to find the
     * corresponding point in material space X and return it.
     *
     * @param c The world space position corresponding to x(X*), where X* is the solution to find.
     * @param Xk The initial guess Xk.
     * @return The material space point X* corresponding to c = x(X*). In other words, return
     * x^{-1}(X).
     */
    Eigen::Vector3d
    inverse(Eigen::Vector3d const& c, Eigen::Vector3d const& Xk, int iterations) const
    {
        int const K              = iterations;
        Eigen::Vector3d const X0 = Xk;
        for (int k = 0; k < K; ++k)
        {
            auto const [gradf, Hf] = hessian_f(c, Xk);
#ifdef _DEBUG
            scalar_type const det = Hf.determinant();
            assert(std::abs(det) > sbs::eps());
#endif
            Eigen::Matrix3d const Hfinv = Hf.inverse();
            Xk                          = Xk - Hfinv * gradf;
        }

#ifdef _DEBUG
        Eigen::Vector3d const x = eval(Xk);
        scalar_type const r     = (x - c).norm();
#endif
        return Xk_plus_1;
    }

    /**
     * @brief
     *
     * @param c
     * @param X
     * @return
     */
    std::pair<Eigen::Vector3d, Eigen::Matrix3d>
    hessian_f(Eigen::Vector3d const& c, Eigen::Vector3d const& X) const
    {
        Eigen::Vector3d gradf{};
        gradf.setZero();
        Eigen::Matrix3d Hf{};
        Hf.setZero();

        auto const [x, gradx, Hx] = hessian_x(X);

        for (int j = 0; j < 3; ++j)
        {
            // gradf contributions
            for (k = 0; k < 3; ++k)
            {
                gradf(j) += (x(k) - c(k)) * gradx[k](j);
            }

            // Hf contributions
            for (int i = 0; i < 3; ++i)
            {
                for (k = 0; k < 3; ++k)
                {
                    scalar_type const c1 = gradx[k](i) * gradx[k](j);
                    scalar_type const c2 = x(k) * Hx(i, j);
                    scalar_type const c3 = c(k) * Hx(i, j);
                    Hf(i, j) += c1 + c2 - c3;
                }
            }
        }

        return std::make_pair(gradf, Hf);
    }

    /**
     * @brief
     * Computes the hessian of the 0th-order reproducible version
     * of this interpolation field.
     *
     * @param X Material space point at which to compute the hessian.
     * @return Tuple (x(X), [grad x_1(X), grad x_2(X), grad x_3(X)], [Hx_1(X), Hx_2(X), Hx_3(X)] )
     */
    std::tuple<Eigen::Vector3d, std::array<Eigen::Vector3d, 3u>, std::array<Eigen::Matrix3d, 3u>>
    hessian_x(Eigen::Vector3d const& X) const
    {
        std::array<Eigen::Vector3d, 3u> gradx{};
        gradx[0].setZero();
        gradx[1].setZero();
        gradx[2].setZero();
        std::array<Eigen::Matrix3d, 3u> Hx{};
        Hx[0].setZero();
        Hx[1].setZero();
        Hx[2].setZero();

        auto const [u, gradu, Hu] = hessian_u(X);
        auto const [v, gradv, Hv] = hessian_v(X);

        for (int k = 0; k < 3; ++k)
        {
            // grad contributions
            scalar_type const v2 = v * v;
            for (int j = 0; j < 3; ++j)
            {
                scalar_type const lhs = gradu[k](j) * v;
                scalar_type const rhs = u(k) * gradv(j);
                gradx[k] += (lhs - rhs) / v2
            }

            // Hessian contributions
            scalar_type const v4 = v2 * v2;
            for (int j = 0; j < 3; ++j)
            {
                for (int i = 0; i < 3; ++i)
                {
                    scalar_type const c1 = (Hu(i, j) * v + gradu[k](j) * gradv(i)) * v2;
                    scalar_type const c2 = (gradu[k](i) * gradv(j) + u(k) * Hv(i, j)) * v2;
                    scalar_type const c3 = (gradu[k](j) * v) * 2. * v * gradv(i);
                    scalar_type const c4 = (u(k) * gradv(j)) * 2. * v * gradv(i);
                    Hx(i, j) += (c1 - c2 - c3 + c4) / v4;
                }
            }
        }

        Eigen::Vector3d const x = u / v;
        return std::make_tuple(x, gradx, Hx);
    }

    /**
     * @brief
     * Returns the hessians of each component of u(X).
     * Also returns the gradients of each component of u(X).
     *
     * @param X Material space position X
     * @return Tuple (u(X), (grad u_1, grad u_2, grad u_3), (Hu_1(X), Hu_2(X), Hu_3(X)))
     */
    std::tuple<Eigen::Vector3d, std::array<Eigen::Vector3d, 3u>, std::array<Eigen::Matrix3d, 3u>>
    hessian_u(Eigen::Vector3d const& X) const
    {
        std::array<Eigen::Vector3d, 3u> gradu{};
        gradu[0].setZero();
        gradu[1].setZero();
        gradu[2].setZero();
        std::array<Eigen::Matrix3d, 3u> Hu{};
        Hu[0].setZero();
        Hu[1].setZero();
        Hu[2].setZero();

        // Compute FEM basis function contribution to gradu and Hu.
        // Linear FEM basis functions do not have a second derivative, thus
        // Hu = 0, and we only compute gradu contributions.
        auto const& cell = (*cells)[e];
        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            auto const i       = cell.node(r);
            bool const has_phi = (*has_basis_function)[i];
            if (has_phi)
            {
                Eigen::Vector3d const& xi     = (*xis)[i];
                auto const& phi               = cell.phi(r);
                Eigen::Vector3d const gradphi = phi.grad(X);

                for (int k = 0; k < 3; ++k)
                {
                    gradu[k] += xi(k) * gradphi;
                    u += xi * phi.eval(X);
                }
            }
        }

        for (index_type const j : js)
        {
            scalar_type const& Vj          = (*Vjs)[j];
            kernel_function_type const& Wj = (*Wjs)[j];
            scalar_type const W            = static_cast<scalar_type>(Wj(X));
            Eigen::Matrix3d const& Fj      = (*Fjs)[j];
            Eigen::Vector3d const& Xj      = (*Xjs)[j];
            Eigen::Vector3d const& xj      = (*xjs)[j];
            Eigen::Vector3d const X_j      = X - Xj;

            Eigen::Vector3d const& gradWj     = Wj.grad(X);
            Eigen::Matrix3d const& hessian_Wj = Wj.hessian(X);

            // Description of indices:
            // j: meshless node indices
            // k: index of component of the vector-valued function u(X)
            // s: index of the component of independent variable X to differentiate against for
            // first-order r: index of the component of independent variable X to differentiate
            // against for second-order m: index of column of deformation gradient in the matrix
            // vector product Fj * (X - Xj)
            Eigen::Matrix3d const dX_dX = Eigen::Matrix3d::Identity();
            for (int k = 0; k < 3; ++k)
            {
                for (int s = 0; s < 3; ++s)
                {
                    // gradu contributions
                    gradu[k](s) += xj(k) * Vj * gradWj(s);
                    for (int m = 0; m < 3; ++m)
                    {
                        gradu[k] += Fj(k, m) * dX_dX(m, s) * Vj * W;
                        gradu[k] += Fj(k, m) * X(m) * Vj * gradWj(s);
                        gradu[k] -= Fj(k, m) * Xj(m) * Vj * gradWj(s);
                    }

                    // Hu contributions
                    for (int r = 0; r < 3; ++r)
                    {
                        Hu[k](r, s) += xj(k) * Vj * hessian_Wj(r, s);
                        for (int m = 0; m < 3; ++m)
                        {
                            Hu[k](r, s) += Fj(k, m) * dX_dX(m, s) * Vj * gradWj(r);
                            Hu[k](r, s) += Fj(k, m) * X(m) * Vj * hessian_Wj(r, s);
                            Hu[k](r, s) -= Fj(k, m) * Xj(m) * Vj * hessian_Wj(r, s);
                        }
                    }
                }
            }
        }

        Eigen::Vector3d const u = get_u(X);
        return std::make_tuple(u, gradu, Hu);
    }

    /**
     * @brief
     * Returns the hessian of v(X).
     * Also returns the gradient of v(X).
     *
     * @param X Material space position X.
     * @return Tuple (v, gradv, Hv)
     */
    std::tuple<scalar_type, Eigen::Vector3d, Eigen::Matrix3d>
    hessian_v(Eigen::Vector3d const& X) const
    {
        Eigen::Vector3d gradv{};
        gradv.setZero();
        Eigen::Matrix3d Hv{};
        Hv.setZero();

        auto const& cell = (*cells)[e];
        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            auto const i       = cell.node(r);
            bool const has_phi = (*has_basis_function)[i];
            if (has_phi)
            {
                auto const& phi               = cell.phi(r);
                Eigen::Vector3d const gradphi = phi.grad(X);
                // gradv contributions
                gradv += gradphi;
            }
        }

        for (index_type const j : js)
        {
            scalar_type const& Vj          = (*Vjs)[j];
            kernel_function_type const& Wj = (*Wjs)[j];
            scalar_type const W            = static_cast<scalar_type>(Wj(X));
            Eigen::Vector3d const& Xj      = (*Xjs)[j];
            Eigen::Vector3d const X_j      = X - Xj;

            Eigen::Vector3d const& gradWj     = Wj.grad(X);
            Eigen::Matrix3d const& hessian_Wj = Wj.hessian(X);

            // gradv contributions
            gradv += Vj * gradWj;

            // Hv contributions
            for (int m = 0; m < 3; ++m)
            {
                for (int n = 0; n < 3; ++n)
                {
                    Hv(m, n) += Vj * hessian_Wj(m, n);
                }
            }
        }

        scalar_type const v = get_v(X);
        return std::make_tuple(v, gradv, Hv);
    }

    Eigen::Vector3d get_u(Eigen::Vector3d const& X) const
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
            xk += (Fj * Xkj + xj) * Vj * Wkj;
        }
        return xk;
    }

    std::vector<Eigen::Matrix3d>* Fjs; ///< Deformation gradients at the meshless nodes
};

template <class FemCellType, class KernelFunctionType = poly6_kernel_t>
struct fem_sph_nodal_deformation_gradient_op_t
{
    using cell_type             = FemCellType;
    using kernel_function_type  = KernelFunctionType;
    using interpolation_op_type = fem_sph_nodal_interpolation_t<cell_type, kernel_function_type>;

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

template <class FemCellType, class KernelFunctionType = poly6_kernel_t>
struct fem_sph_efg_deformation_gradient_op_t
{
    using cell_type             = FemCellType;
    using kernel_function_type  = KernelFunctionType;
    using interpolation_op_type = fem_sph_efg_interpolation_t<cell_type, kernel_function_type>;

    fem_sph_efg_deformation_gradient_op_t(
        std::vector<Eigen::Vector3d> const* Xks,
        index_type k,
        interpolation_op_type const& interpolation_op)
        : Xks(Xks),
          k(k),
          interpolation_op(interpolation_op),
          phi_is(),
          phi_js(),
          gradphi_is(),
          gradphi_js(),
          corrected_gradphi_is(),
          corrected_gradphi_js(),
          Lk()
    {
        // precompute correction matrix L and shape function derivatives
        store_L_and_grads();
    }

    void store_L_and_grads()
    {
        phi_is.clear();
        phi_js.clear();
        gradphi_is.clear();
        gradphi_js.clear();
        corrected_gradphi_is.clear();
        corrected_gradphi_js.clear();
        Lk.setZero();

        auto const& Xis                = *interpolation_op.Xis;
        auto const& cells              = *interpolation_op.cells;
        auto const& cell               = cells[interpolation_op.e];
        auto const& has_basis_function = *interpolation_op.has_basis_function;

        Eigen::Vector3d const& Xk = (*Xks)[k];

        // Precompute gradients
        scalar_type sum_phi_i_phi_j = 0.;
        Eigen::Vector3d sum_gradphi_i_gradphi_j{0., 0., 0.};

        for (auto r = 0u; r < cell.node_count(); ++r)
        {
            auto const i       = cell.node(r);
            bool const has_phi = has_basis_function[i];
            if (has_phi)
            {
                auto const& phi               = cell.phi(r);
                Eigen::Vector3d const gradphi = phi.grad(Xk);

                Eigen::Vector3d const& Xi = Xis[i];
                scalar_type const phi_i   = phi.eval(Xk);

                phi_is.push_back(phi_i);
                gradphi_is.push_back(gradphi);
                sum_phi_i_phi_j += phi_i;
                sum_gradphi_i_gradphi_j += gradphi;
            }
        }

        std::vector<scalar_type> const& Vjs          = *interpolation_op.Vjs;
        std::vector<kernel_function_type> const& Wjs = *interpolation_op.Wjs;
        std::vector<Eigen::Vector3d> const& Xjs      = *interpolation_op.Xjs;

        auto const& js = interpolation_op.js;
        for (index_type const j : js)
        {
            scalar_type const& Vj         = Vjs[j];
            kernel_function_type const& W = Wjs[j];
            scalar_type const Wj          = W.eval(Xk);

            scalar_type const phi_j       = Vj * Wj;
            Eigen::Vector3d const gradWkj = W.grad(Xk);
            Eigen::Vector3d const gradphi = Vj * gradWkj;

            phi_js.push_back(phi_j);
            gradphi_js.push_back(gradphi);
            sum_phi_i_phi_j += phi_j;
            sum_gradphi_i_gradphi_j += gradphi;
        }

        // Correct the gradients
        scalar_type const squared_sum_phi_i_phi_j = (sum_phi_i_phi_j * sum_phi_i_phi_j);
        scalar_type const shepard_coefficient     = 1. / squared_sum_phi_i_phi_j;

        for (auto a = 0u; a < gradphi_is.size(); ++a)
        {
            Eigen::Vector3d const corrected_gradphi_i =
                shepard_coefficient *
                ((gradphi_is[a] * sum_phi_i_phi_j) - (phi_is[a] * sum_gradphi_i_gradphi_j));

            corrected_gradphi_is.push_back(corrected_gradphi_i);
        }
        for (auto a = 0u; a < gradphi_js.size(); ++a)
        {
            Eigen::Vector3d const corrected_gradphi_j =
                shepard_coefficient *
                ((gradphi_js[a] * sum_phi_i_phi_j) - (phi_js[a] * sum_gradphi_i_gradphi_j));

            corrected_gradphi_js.push_back(corrected_gradphi_j);
        }

        // Compute correction matrix L
        Eigen::Matrix3d Linv{};
        Linv.setZero();
        assert(corrected_gradphi_is.size() == interpolation_op.is.size());
        for (auto a = 0u; a < gradphi_is.size(); ++a)
        {
            index_type const i                       = interpolation_op.is[a];
            Eigen::Vector3d const& Xi                = Xis[i];
            Eigen::Vector3d const& corrected_gradphi = corrected_gradphi_is[a];
            Eigen::Vector3d const Xik                = Xi - Xk;
            Linv += corrected_gradphi * Xik.transpose();
        }
        assert(corrected_gradphi_js.size() == interpolation_op.js.size());
        for (auto a = 0u; a < corrected_gradphi_js.size(); ++a)
        {
            index_type const j                       = js[a];
            Eigen::Vector3d const& Xj                = Xjs[j];
            Eigen::Vector3d const& corrected_gradphi = corrected_gradphi_js[a];
            Eigen::Vector3d const Xjk                = Xj - Xk;
            Linv += corrected_gradphi * Xjk.transpose();
        }
#ifdef _DEBUG
        scalar_type const det = Linv.determinant();
        assert(std::abs(det) > sbs::eps());
#endif
        Lk = Linv.inverse();
    }

    Eigen::Matrix3d eval() const
    {
        std::vector<index_type> const& is       = interpolation_op.is;
        std::vector<index_type> const& js       = interpolation_op.js;
        std::vector<Eigen::Vector3d> const& xis = *interpolation_op.xis;
        std::vector<Eigen::Vector3d> const& xjs = *interpolation_op.xjs;

        Eigen::Matrix3d F{};
        F.setZero();

        for (auto a = 0u; a < is.size(); ++a)
        {
            index_type const i                       = is[a];
            Eigen::Vector3d const& corrected_gradphi = corrected_gradphi_is[a];
            Eigen::Vector3d const& xi                = xis[i];
            F += xi * (Lk * corrected_gradphi).transpose();
        }
        for (auto a = 0u; a < js.size(); ++a)
        {
            index_type const j                       = js[a];
            Eigen::Vector3d const& corrected_gradphi = corrected_gradphi_js[a];
            Eigen::Vector3d const& xj                = xjs[j];
            F += xj * (Lk * corrected_gradphi).transpose();
        }

        return F;
    }

    /**
     * @brief
     * Returns a pair of the gradients w.r.t. fem dofs and gradients w.r.t. to sph dofs.
     * @return ([grad \phi_is], [grad \phi_js])
     */
    std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> dFdx() const
    {
        std::vector<index_type> const& is = interpolation_op.is;
        std::vector<index_type> const& js = interpolation_op.js;

        std::vector<Eigen::Vector3d> dFdxis{};
        dFdxis.resize(is.size(), Eigen::Vector3d{0., 0., 0.});
        std::vector<Eigen::Vector3d> dFdxjs{};
        dFdxjs.resize(js.size(), Eigen::Vector3d{0., 0., 0.});

        for (auto a = 0u; a < is.size(); ++a)
        {
            index_type const i                       = is[a];
            Eigen::Vector3d const& corrected_gradphi = corrected_gradphi_is[a];
            dFdxis[a]                                = Lk * corrected_gradphi;
        }
        for (auto a = 0u; a < js.size(); ++a)
        {
            index_type const j                      = js[a];
            Eigen::Vector3d const corrected_gradphi = corrected_gradphi_js[a];

            dFdxjs[a] = Lk * corrected_gradphi;
        }

        return std::make_pair(dFdxis, dFdxjs);
    }

    std::vector<Eigen::Vector3d> const* Xks; ///< EFG integration points

    index_type k;                                      ///< Index of the sph node
    interpolation_op_type const& interpolation_op;     ///< The mixed sph+fem interpolation
    std::vector<scalar_type> phi_is;                   ///< Precomputed FEM basis functions
    std::vector<scalar_type> phi_js;                   ///< Precomputed SPH basis functions
    std::vector<Eigen::Vector3d> gradphi_is;           ///< Cached fem basis function gradients
                                                       ///< in this interpolation
    std::vector<Eigen::Vector3d> gradphi_js;           ///< Cached sph basis function gradients
                                                       ///< in this interpolation
    std::vector<Eigen::Vector3d> corrected_gradphi_is; ///< Corrected FEM basis function gradients
    std::vector<Eigen::Vector3d> corrected_gradphi_js; ///< Corrected SPH basis function gradients
    Eigen::Matrix3d Lk;                                ///< Lsph^-1 * Lfem
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_FEM_SPH_H
