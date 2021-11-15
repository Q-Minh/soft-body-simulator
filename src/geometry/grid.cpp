#include "sbs/geometry/grid.h"

#include "sbs/aliases.h"

#include <numeric>

namespace sbs {
namespace geometry {

grid_t::grid_t(tetrahedral_domain_t const& domain, Eigen::Vector3i const& dims)
    : bounding_box_(), dims_(dims)
{
    auto constexpr min  = std::numeric_limits<scalar_type>::lowest();
    auto constexpr max  = std::numeric_limits<scalar_type>::max();
    bounding_box_.min() = Eigen::Vector3d{min, min, min};
    bounding_box_.max() = Eigen::Vector3d{max, max, max};

    auto const N = domain.topology().vertex_count();
    for (auto i = 0u; i < N; ++i)
    {
        Eigen::Vector3d const& Xi = domain.position(i);
        bounding_box_.min()       = bounding_box_.min().cwiseMin(Xi);
        bounding_box_.max()       = bounding_box_.max().cwiseMax(Xi);
    }
}

void grid_t::walk(std::function<void(Eigen::Vector3d const&)> const& f) const
{
    Eigen::Vector3d const dX      = delta();
    Eigen::Vector3d const half_dX = dX / 2.;

    for (int k = 0; k < dims_(2); ++k)
    {
        for (int j = 0; j < dims_(1); ++j)
        {
            for (int i = 0; i < dims_(0); ++i)
            {
                Eigen::Vector3i const ijk{i, j, k};
                Eigen::Vector3d const grid_node =
                    bounding_box_.min() +
                    Eigen::Vector3d{dX.array() * ijk.cast<scalar_type>().array()};

                Eigen::Vector3d const center = grid_node + half_dX;

                f(center);
            }
        }
    }
}

Eigen::Vector3d grid_t::delta() const
{
    Eigen::Vector3d const d =
        (bounding_box_.max() - bounding_box_.min()).array() / dims_.cast<scalar_type>().array();
    return d;
}

} // namespace geometry
} // namespace sbs
