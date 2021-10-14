#include "sbs/physics/xpbd/green_constraint.h"

#include <Eigen/LU>
#include <Eigen/SVD>
#include <sbs/physics/simulation.h>

namespace sbs {
namespace physics {
namespace xpbd {

green_constraint_t::green_constraint_t(
    scalar_type const alpha,
    scalar_type const beta,
    simulation_t const& simulation,
    index_type bi,
    index_type v1,
    index_type v2,
    index_type v3,
    index_type v4,
    scalar_type young_modulus,
    scalar_type poisson_ratio)
    : constraint_t(alpha, beta),
      bi_(bi),
      v1_(v1),
      v2_(v2),
      v3_(v3),
      v4_(v4),
      DmInv_(),
      V0_(),
      mu_(),
      lambda_()
{
    auto const& p1 = simulation.particles()[bi_][v1_];
    auto const& p2 = simulation.particles()[bi_][v2_];
    auto const& p3 = simulation.particles()[bi_][v3_];
    auto const& p4 = simulation.particles()[bi_][v4_];

    Eigen::Matrix3d Dm;
    Dm.col(0) = (p1.x0() - p4.x0()).transpose();
    Dm.col(1) = (p2.x0() - p4.x0()).transpose();
    Dm.col(2) = (p3.x0() - p4.x0()).transpose();

    DmInv_  = Dm.inverse();
    V0_     = (1. / 6.) * Dm.determinant();
    mu_     = (young_modulus) / (2. * (1 + poisson_ratio));
    lambda_ = (young_modulus * poisson_ratio) / ((1 + poisson_ratio) * (1 - 2 * poisson_ratio));
}

void green_constraint_t::project_positions(simulation_t& simulation, scalar_type dt)
{
    auto& p1 = simulation.particles()[bi_][v1_];
    auto& p2 = simulation.particles()[bi_][v2_];
    auto& p3 = simulation.particles()[bi_][v3_];
    auto& p4 = simulation.particles()[bi_][v4_];

    scalar_type const w1 = p1.invmass();
    scalar_type const w2 = p2.invmass();
    scalar_type const w3 = p3.invmass();
    scalar_type const w4 = p4.invmass();

    auto const Vsigned        = signed_volume(p1.xi(), p2.xi(), p3.xi(), p4.xi());
    bool const is_V_positive  = Vsigned >= 0.;
    bool const is_V0_positive = V0_ >= 0.;
    bool const is_tet_inverted =
        (is_V_positive && !is_V0_positive) || (!is_V_positive && is_V0_positive);

    scalar_type constexpr epsilon = 1e-20;

    Eigen::Matrix3d Ds;
    Ds.col(0) = (p1.xi() - p4.xi());
    Ds.col(1) = (p2.xi() - p4.xi());
    Ds.col(2) = (p3.xi() - p4.xi());

    Eigen::Matrix3d const F = Ds * DmInv_;
    Eigen::Matrix3d const I = Eigen::Matrix3d::Identity();

    // TODO: Implement correct inversion handling described in
    // Irving, Geoffrey, Joseph Teran, and Ronald Fedkiw. "Invertible finite elements for robust
    // simulation of large deformation." Proceedings of the 2004 ACM SIGGRAPH/Eurographics symposium
    // on Computer animation. 2004.
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
    scalar_type const dt2         = dt * dt;
    scalar_type const alpha_tilde = alpha() / dt2;
    scalar_type const beta_tilde  = beta() * dt2;
    scalar_type const gamma       = alpha_tilde * beta_tilde / dt;

    // clang-format off
    scalar_type const gradC_dot_displacement =
        f1.dot(p1.xi() - p1.xn()) + 
        f2.dot(p2.xi() - p2.xn()) + 
        f3.dot(p3.xi() - p3.xn()) +
        f4.dot(p4.xi() - p4.xn());
    // clang-format on

    scalar_type const delta_lagrange_num =
        -(C + alpha_tilde * lagrange_) + gamma * gradC_dot_displacement;
    scalar_type const delta_lagrange_den = (1. + gamma) * (weighted_sum_of_gradients) + alpha_tilde;
    scalar_type const delta_lagrange     = delta_lagrange_num / delta_lagrange_den;

    lagrange_ += delta_lagrange;
    // because f = - grad(potential), then grad(potential) = -f and thus grad(C) = -f
    p1.xi() += w1 * -f1 * delta_lagrange;
    p2.xi() += w2 * -f2 * delta_lagrange;
    p3.xi() += w3 * -f3 * delta_lagrange;
    p4.xi() += w4 * -f4 * delta_lagrange;
}

scalar_type green_constraint_t::signed_volume(
    Eigen::Vector3d const& p1,
    Eigen::Vector3d const& p2,
    Eigen::Vector3d const& p3,
    Eigen::Vector3d const& p4) const
{
    Eigen::Matrix3d Dm;
    Dm.col(0) = (p1 - p4);
    Dm.col(1) = (p2 - p4);
    Dm.col(2) = (p3 - p4);

    scalar_type const V = (1. / 6.) * Dm.determinant();
    return V;
}

} // namespace xpbd
} // namespace physics
} // namespace sbs