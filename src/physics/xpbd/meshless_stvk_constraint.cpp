#include <numeric>
#include <sbs/physics/simulation.h>
#include <sbs/physics/xpbd/meshless_stvk_constraint.h>

namespace sbs {
namespace physics {
namespace xpbd {

meshless_stvk_constraint_t::meshless_stvk_constraint_t(
    scalar_type const alpha,
    scalar_type const beta,
    simulation_t const& simulation,
    index_type bi,
    index_type vi,
    scalar_type young_modulus,
    scalar_type poisson_ratio,
    mechanics::meshless_node_t const& node)
    : constraint_t(alpha, beta), node_(node), bi_(bi), vi_(vi), mu_(), lambda_()
{
    mu_     = (young_modulus) / (2. * (1 + poisson_ratio));
    lambda_ = (young_modulus * poisson_ratio) / ((1 + poisson_ratio) * (1 - 2 * poisson_ratio));
}

void meshless_stvk_constraint_t::project_positions(simulation_t& simulation, scalar_type dt)
{
    auto& particles                           = simulation.particles()[bi_];
    std::vector<index_type> const& neighbours = node_.neighbours();
    functions::poly6_kernel_t const& kernel   = node_.kernel();
    Eigen::Matrix3d Fi{};
    Fi.setZero();
    for (std::size_t k = 0u; k < neighbours.size(); ++k)
    {
        index_type const j             = neighbours[k];
        scalar_type const Vj           = node_.Vjs()[k];
        Eigen::Matrix3d const& Li      = node_.Li();
        Eigen::Vector3d const& gradWij = node_.gradWij()[k];
        particle_t const& pi           = particles[vi_];
        particle_t const& pj           = particles[j];
        Eigen::Vector3d const xji      = pj.xi() - pi.xi();
        Eigen::Matrix3d const Fij      = Vj * xji * (Li * gradWij).transpose();
        Fi += Fij;
    }

    Eigen::Matrix3d const Ei =
        scalar_type{0.5} * (Fi.transpose() * Fi) - Eigen::Matrix3d::Identity();
    scalar_type const Eitrace = Ei.trace();
    scalar_type const Psi =
        mu_ * (Ei.array() * Ei.array()).sum() + 0.5 * lambda_ * (Eitrace * Eitrace);
    scalar_type const C = node_.Vi() * Psi;

    // d (Psi) / d (Fi) yields a 3x3 matrix
    Eigen::Matrix3d const dPsidFi =
        Fi * (2. * mu_ * Ei + lambda_ * Eitrace * Eigen::Matrix3d::Identity());

    // d (Fi) / d (xk) yields a 3x3x3 tensor, but each d (Fi) / d (xk)^(l), for l=1,2,3,
    // l representing the component of xk, d (Fi) / d (xk)^(l) is a 3x3 matrix
    std::vector<Eigen::Vector3d> gradC{};
    gradC.resize(neighbours.size(), Eigen::Vector3d{0., 0., 0.});
    index_type const i_idx =
        std::distance(neighbours.begin(), std::find(neighbours.begin(), neighbours.end(), vi_));
    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        index_type const k = neighbours[a];

        if (k == vi_)
            continue;

        Eigen::Matrix3d const& Li = node_.Li();

        index_type const j              = k;
        scalar_type const Vj            = node_.Vjs()[a];
        Eigen::Vector3d const& gradWij  = node_.gradWij()[a];
        Eigen::Vector3d const LiGradWij = Li * gradWij;

        Eigen::Vector3d const& dFidxj_1 = Vj * LiGradWij;
        Eigen::Vector3d const& dFidxj_2 = dFidxj_1;
        Eigen::Vector3d const& dFidxj_3 = dFidxj_1;

        Eigen::Vector3d gradPsi{0., 0., 0.};
        gradPsi(0u) = node_.Vi() * dPsidFi.row(0u).dot(dFidxj_1);
        gradPsi(1u) = node_.Vi() * dPsidFi.row(1u).dot(dFidxj_2);
        gradPsi(2u) = node_.Vi() * dPsidFi.row(2u).dot(dFidxj_3);

        gradC[a] += gradPsi;
        gradC[i_idx] -= gradPsi;
    }

    scalar_type weighted_sum_of_gradients{0.};
    scalar_type gradC_dot_displacement{0.};
    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        index_type const j   = neighbours[a];
        particle_t const& pj = particles[j];
        weighted_sum_of_gradients += pj.invmass() * gradC[a].squaredNorm();
        gradC_dot_displacement += gradC[a].dot(pj.xi() - pj.xn());
    }

    scalar_type constexpr epsilon = 1e-20;
    if (weighted_sum_of_gradients < epsilon)
        return;

    scalar_type const dt2         = dt * dt;
    scalar_type const alpha_tilde = alpha() / dt2;
    scalar_type const beta_tilde  = beta() * dt2;
    scalar_type const gamma       = alpha_tilde * beta_tilde / dt;

    scalar_type const delta_lagrange_num =
        -(C + alpha_tilde * lagrange_) + gamma * gradC_dot_displacement;
    scalar_type const delta_lagrange_den = (1. + gamma) * (weighted_sum_of_gradients) + alpha_tilde;
    scalar_type const delta_lagrange     = delta_lagrange_num / delta_lagrange_den;

    lagrange_ += delta_lagrange;

    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        index_type const j = neighbours[a];
        particle_t& pj     = particles[j];
        pj.xi() += pj.invmass() * gradC[a] * delta_lagrange;
    }
}

} // namespace xpbd
} // namespace physics
} // namespace sbs