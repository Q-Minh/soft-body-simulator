#include <Eigen/LU>
#include <cassert>
#include <iterator>
#include <sbs/physics/functions/kernel.h>
#include <sbs/physics/mechanics/meshless_sph_node.h>
#include <sbs/physics/particle.h>

namespace sbs {
namespace physics {
namespace mechanics {

meshless_sph_node_t::meshless_sph_node_t(
    index_type const i,
    functions::poly6_kernel_t const& kernel)
    : ni_(i), neighbours_(), Wij_(), gradWij_(), Vjs_(), Vi_(0.), Ci_(), kernel_(kernel), Fi_(), xi_()
{
    Fi_.setIdentity();
    xi_.setZero();
}

void meshless_sph_node_t::initialize(
    Eigen::Vector3d const& xi,
    std::vector<Eigen::Vector3d const*> const& pj,
    std::vector<index_type> const& neighbours)
{
    assert(pj.size() == neighbours.size());

    neighbours_ = neighbours;

    functions::poly6_kernel_t& W = kernel_;

    for (std::size_t k = 0u; k < neighbours.size(); ++k)
    {
        Eigen::Vector3d const Xj = *pj[k];
        // W(||Xi - Xj||)
        scalar_type const Wi = W(Xj);
        // grad W(||Xi - Xj||)
        Eigen::Vector3d const gradW = W.grad(Xj);

        Wij_.push_back(Wi);
        gradWij_.push_back(gradW);

        Vi_ += Wi;
    }

    // Vi = 1 / (sum W(||Xi - Xj||))
    Vi_ = 1. / Vi_;

    Fi_.setIdentity();
    xi_ = xi;
}

index_type meshless_sph_node_t::Ni() const
{
    return ni_;
}

std::vector<index_type> const& meshless_sph_node_t::neighbours() const
{
    return neighbours_;
}

std::vector<scalar_type> const& meshless_sph_node_t::Wij() const
{
    return Wij_;
}

std::vector<Eigen::Vector3d> const& meshless_sph_node_t::gradWij() const
{
    return gradWij_;
}

std::vector<scalar_type> const& meshless_sph_node_t::Vjs() const
{
    return Vjs_;
}

scalar_type meshless_sph_node_t::Vi() const
{
    return Vi_;
}

Eigen::Matrix3d const& meshless_sph_node_t::Li() const
{
    return Ci_;
}

functions::poly6_kernel_t const& meshless_sph_node_t::kernel() const
{
    return kernel_;
}

Eigen::Vector3d const& meshless_sph_node_t::Xi() const
{
    return kernel_.xi();
}

Eigen::Vector3d& meshless_sph_node_t::Xi()
{
    return kernel_.xi();
}

Eigen::Matrix3d const& meshless_sph_node_t::Fi() const
{
    return Fi_;
}

Eigen::Matrix3d& meshless_sph_node_t::Fi()
{
    return Fi_;
}

Eigen::Vector3d const& meshless_sph_node_t::xi() const
{
    return xi_;
}

Eigen::Vector3d& meshless_sph_node_t::xi()
{
    return xi_;
}

Eigen::Matrix3d meshless_sph_node_t::Ei() const
{
    Eigen::Matrix3d const Ei =
        scalar_type{0.5} * ((Fi_.transpose() * Fi_) - Eigen::Matrix3d::Identity());
    return Ei;
}

std::pair<scalar_type, Eigen::Matrix3d> meshless_sph_node_t::Psi_dPsidFi(
    Eigen::Matrix3d const& Ei,
    scalar_type const mu,
    scalar_type const lambda) const
{
    scalar_type const Eitrace = Ei.trace();
    scalar_type const Psi =
        mu * (Ei.array() * Ei.array()).sum() + 0.5 * lambda * (Eitrace * Eitrace);
    scalar_type const C = Vi_ * Psi;

    // d (Psi) / d (Fi) yields a 3x3 matrix
    Eigen::Matrix3d const dPsidFi =
        Fi_ * (2. * mu * Ei + lambda * Eitrace * Eigen::Matrix3d::Identity());

    return std::make_pair(Psi, dPsidFi);
}

scalar_type const meshless_sph_node_t::E(scalar_type const Psi) const
{
    return Vi_ * Psi;
}

std::vector<Eigen::Vector3d> meshless_sph_node_t::dEdxk(Eigen::Matrix3d const& dPsidFi) const
{
    std::vector<index_type> const& neighbours = neighbours_;

    std::vector<Eigen::Vector3d> gradC{};
    gradC.resize(neighbours.size(), Eigen::Vector3d{0., 0., 0.});
    index_type const i_idx = static_cast<index_type>(
        std::distance(neighbours_.begin(), std::find(neighbours_.begin(), neighbours_.end(), ni_)));
    for (std::size_t a = 0u; a < neighbours.size(); ++a)
    {
        index_type const k = neighbours[a];

        if (k == ni_)
            continue;

        Eigen::Matrix3d const& Li = Ci_;

        index_type const j              = k;
        scalar_type const Vj            = Vjs_[a];
        Eigen::Vector3d const& gradWij  = gradWij_[a];
        Eigen::Vector3d const LiGradWij = Li * gradWij;

        Eigen::Vector3d const& dFidxj_1 = Vj * LiGradWij;
        Eigen::Vector3d const& dFidxj_2 = dFidxj_1;
        Eigen::Vector3d const& dFidxj_3 = dFidxj_1;

        Eigen::Vector3d gradPsi{0., 0., 0.};
        gradPsi(0u) = Vi_ * dPsidFi.row(0u).dot(dFidxj_1);
        gradPsi(1u) = Vi_ * dPsidFi.row(1u).dot(dFidxj_2);
        gradPsi(2u) = Vi_ * dPsidFi.row(2u).dot(dFidxj_3);

        gradC[a] += gradPsi;
        gradC[i_idx] -= gradPsi;
    }

    return gradC;
}

void meshless_sph_node_t::cache_Li_Vj(std::vector<meshless_sph_node_t const*> const& neighbours)
{
    Eigen::Matrix3d Li{};
    Li.setZero();

    for (std::size_t j = 0u; j < neighbours.size(); ++j)
    {
        Eigen::Vector3d const Xji     = neighbours[j]->Xi() - Xi();
        scalar_type const Vj          = neighbours[j]->Vi();
        Eigen::Vector3d const gradWij = gradWij_[j];
        Eigen::Matrix3d const L       = Vj * gradWij * Xji.transpose();
        Li += L;

        Vjs_.push_back(Vj);
    }

    Eigen::Matrix3d LiInv{};
    bool invertible{false};
    double constexpr error = 1e-18;
    Li.computeInverseWithCheck(LiInv, invertible, error);
    assert(invertible);
    Ci_ = LiInv;
}

} // namespace mechanics
} // namespace physics
} // namespace sbs