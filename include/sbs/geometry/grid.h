#ifndef SBS_GEOMETRY_GRID_H
#define SBS_GEOMETRY_GRID_H

#include "tetrahedral_domain.h"

#include <Eigen/Geometry>
#include <functional>

namespace sbs {
namespace geometry {

class grid_t
{
  public:
    grid_t() = default;
    grid_t(Eigen::AlignedBox3d const& aabb, Eigen::Vector3i const& dims)
        : bounding_box_(aabb), dims_(dims)
    {
    }
    grid_t(tetrahedral_domain_t const& domain, Eigen::Vector3i const& dims);

    /**
     * @brief Walk over grid cell centers
     * @param f Handler to call on every grid cell center
     */
    void walk(std::function<void(Eigen::Vector3d const&)> const& f) const;

    Eigen::Vector3d delta() const;

  private:
    Eigen::AlignedBox3d bounding_box_;
    Eigen::Vector3i dims_;
};

} // namespace geometry
} // namespace sbs

#endif // SBS_GEOMETRY_GRID_H