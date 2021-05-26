#include "physics/xpbd/green_constraint.h"

#include <Eigen/SVD>
#include <Eigen/LU>

namespace sbs {
namespace physics {
namespace xpbd {

green_constraint_t::green_constraint_t(
    scalar_type const alpha,
    positions_type const& positions,
    index_pair_type const& vb1,
    index_pair_type const& vb2,
    index_pair_type const& vb3,
    index_pair_type const& vb4,
    scalar_type young_modulus,
    scalar_type poisson_ratio)
    : constraint_t(alpha),
      v1_(vb1.first),
      v2_(vb2.first),
      v3_(vb3.first),
      v4_(vb4.first),
      b1_(vb1.second),
      b2_(vb2.second),
      b3_(vb3.second),
      b4_(vb4.second),
      DmInv_(),
      V0_(),
      mu_(),
      lambda_()
{
    Eigen::Vector3d const p1 = positions.col(v1_);
    Eigen::Vector3d const p2 = positions.col(v2_);
    Eigen::Vector3d const p3 = positions.col(v3_);
    Eigen::Vector3d const p4 = positions.col(v4_);

    Eigen::Matrix3d Dm;
    Dm.col(0) = (p1 - p4).transpose();
    Dm.col(1) = (p2 - p4).transpose();
    Dm.col(2) = (p3 - p4).transpose();

    DmInv_  = Dm.inverse();
    V0_     = (1. / 6.) * Dm.determinant();
    mu_     = (young_modulus) / (2. * (1 + poisson_ratio));
    lambda_ = (young_modulus * poisson_ratio) / ((1 + poisson_ratio) * (1 - 2 * poisson_ratio));
}

void green_constraint_t::project(
    std::vector<positions_type>& positions,
    std::vector<masses_type> const& masses,
    scalar_type& lagrange_multiplier,
    scalar_type const dt) const
{
    auto const p1 = positions[b1_].col(v1_);
    auto const p2 = positions[b2_].col(v2_);
    auto const p3 = positions[b3_].col(v3_);
    auto const p4 = positions[b4_].col(v4_);

    auto const w1 = 1. / masses[b1_](v1_);
    auto const w2 = 1. / masses[b2_](v2_);
    auto const w3 = 1. / masses[b3_](v3_);
    auto const w4 = 1. / masses[b4_](v4_);

    auto const Vsigned        = signed_volume(p1, p2, p3, p4);
    bool const is_V_positive  = Vsigned >= 0.;
    bool const is_V0_positive = V0_ >= 0.;
    bool const is_tet_inverted =
        (is_V_positive && !is_V0_positive) || (!is_V_positive && is_V0_positive);
    scalar_type constexpr epsilon = 1e-20;

    Eigen::Matrix3d Ds;
    Ds.col(0) = (p1 - p4).transpose();
    Ds.col(1) = (p2 - p4).transpose();
    Ds.col(2) = (p3 - p4).transpose();

    Eigen::Matrix3d const F = Ds * DmInv_;
    Eigen::Matrix3d const I = Eigen::Matrix3d::Identity();

    // TODO: Implement correct inversion handling described in
    // Irving, Geoffrey, Joseph Teran, and Ronald Fedkiw. "Invertible finite elements for robust
    // simulation of large deformation." Proceedings of the 2004 ACM SIGGRAPH/Eurographics symposium
    // on Computer animation. 2004.
    // scalar_type psi{};
    // Eigen::Matrix3d Piola;
    Eigen::JacobiSVD<Eigen::Matrix3d> UFhatV(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d const Fsigma = UFhatV.singularValues();
    Eigen::Matrix3d Fhat;
    Fhat.setZero();
    Fhat(0, 0) = Fsigma(0);
    Fhat(1, 1) = Fsigma(1);
    Fhat(2, 2) = Fsigma(2);

    Eigen::Matrix3d U       = UFhatV.matrixU();
    Eigen::Matrix3d const V = UFhatV.matrixV();

    if (is_tet_inverted)
    {
        Fhat(2, 2) = -Fhat(2, 2);
        U.col(2)   = -U.col(2);
    }

    // stress reaches maximum at 58% compression
    scalar_type constexpr min_singular_value = 0.577;
    Fhat(0, 0)                               = std::max(Fhat(0, 0), min_singular_value);
    Fhat(1, 1)                               = std::max(Fhat(1, 1), min_singular_value);
    Fhat(2, 2)                               = std::max(Fhat(2, 2), min_singular_value);

    Eigen::Matrix3d const Ehat     = 0.5 * (Fhat.transpose() * Fhat - I);
    scalar_type const EhatTrace    = Ehat.trace();
    Eigen::Matrix3d const Piolahat = Fhat * ((2. * mu_ * Ehat) + (lambda_ * EhatTrace * I));

    Eigen::Matrix3d const E  = U * Ehat * V.transpose();
    scalar_type const Etrace = E.trace();
    scalar_type const psi = mu_ * (E.array() * E.array()).sum() + 0.5 * lambda_ * Etrace * Etrace;

    Eigen::Matrix3d const Piola = U * Piolahat * V.transpose();

    // H is the negative gradient of the elastic potential
    scalar_type const V0     = std::abs(V0_);
    Eigen::Matrix3d const H  = -V0 * Piola * DmInv_.transpose();
    Eigen::Vector3d const f1 = H.col(0);
    Eigen::Vector3d const f2 = H.col(1);
    Eigen::Vector3d const f3 = H.col(2);
    Eigen::Vector3d const f4 = -(f1 + f2 + f3);

    // clang-format off
     auto const weighted_sum_of_gradients =
        w1 * f1.squaredNorm() +
        w2 * f2.squaredNorm() +
        w3 * f3.squaredNorm() +
        w4 * f4.squaredNorm();
    // clang-format on

    if (weighted_sum_of_gradients < epsilon)
        return;

    scalar_type const C           = V0 * psi;
    scalar_type const alpha_tilde = alpha_ / (dt * dt);
    scalar_type const delta_lagrange =
        -(C + alpha_tilde * lagrange_multiplier) / (weighted_sum_of_gradients + alpha_tilde);

    lagrange_multiplier += delta_lagrange;
    // because f = - grad(potential), then grad(potential) = -f and thus grad(C) = -f
    positions[b1_].col(v1_) += w1 * -f1 * delta_lagrange;
    positions[b2_].col(v2_) += w2 * -f2 * delta_lagrange;
    positions[b3_].col(v3_) += w3 * -f3 * delta_lagrange;
    positions[b4_].col(v4_) += w4 * -f4 * delta_lagrange;
}

green_constraint_t::scalar_type green_constraint_t::signed_volume(
    Eigen::Vector3d const& p1,
    Eigen::Vector3d const& p2,
    Eigen::Vector3d const& p3,
    Eigen::Vector3d const& p4) const
{
    Eigen::Matrix3d Dm;
    Dm.col(0) = (p1 - p4).transpose();
    Dm.col(1) = (p2 - p4).transpose();
    Dm.col(2) = (p3 - p4).transpose();

    scalar_type const V = (1. / 6.) * Dm.determinant();
    return V;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs