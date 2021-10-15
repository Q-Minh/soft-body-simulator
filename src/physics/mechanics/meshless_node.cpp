#include <Eigen/LU>
#include <cassert>
#include <sbs/physics/functions/kernel.h>
#include <sbs/physics/mechanics/meshless_node.h>
#include <sbs/physics/particle.h>

namespace sbs {
namespace physics {
namespace mechanics {

meshless_node_t::meshless_node_t(index_type const i, functions::poly6_kernel_t const& kernel)
    : ni_(i), neighbours_(), Wij_(), gradWij_(), Vjs_(), Vi_(), Ci_(), kernel_(kernel)
{
}

meshless_node_t::meshless_node_t(
    index_type const i,
    std::vector<Eigen::Vector3d const*> const& pj,
    std::vector<index_type> const& neighbours,
    functions::poly6_kernel_t const& kernel)
    : ni_(i), neighbours_(neighbours), Wij_(), gradWij_(), Vjs_(), Vi_(0.), Ci_(), kernel_(kernel)
{
    assert(pj.size() == neighbours.size());

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
}

index_type meshless_node_t::Ni() const
{
    return ni_;
}

std::vector<index_type> const& meshless_node_t::neighbours() const
{
    return neighbours_;
}

std::vector<scalar_type> const& meshless_node_t::Wij() const
{
    return Wij_;
}

std::vector<Eigen::Vector3d> const& meshless_node_t::gradWij() const
{
    return gradWij_;
}

std::vector<scalar_type> const& meshless_node_t::Vjs() const
{
    return Vjs_;
}

scalar_type meshless_node_t::Vi() const
{
    return Vi_;
}

Eigen::Matrix3d const& meshless_node_t::Li() const
{
    return Ci_;
}

functions::poly6_kernel_t const& meshless_node_t::kernel() const
{
    return kernel_;
}

Eigen::Vector3d const& meshless_node_t::Xi() const
{
    return kernel_.xi();
}

void meshless_node_t::cache_Li_Vj(std::vector<meshless_node_t const*> const& neighbours)
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