#ifndef SBS_GEOMETRY_TRANSFORM_H
#define SBS_GEOMETRY_TRANSFORM_H

#include "sbs/aliases.h"

#include <Eigen/Core>
#include <vector>

namespace sbs {
namespace geometry {

std::vector<Eigen::Vector3d>
translate(std::vector<Eigen::Vector3d> const& points, Eigen::Vector3d const& translation);

Eigen::Vector3d center_of_geometry(std::vector<Eigen::Vector3d> const& points);

std::vector<Eigen::Vector3d> rotate(
    std::vector<Eigen::Vector3d> const& points,
    Eigen::AngleAxis<scalar_type> const& angle_axis,
    Eigen::Vector3d const& origin);

std::vector<Eigen::Vector3d> scale(
    std::vector<Eigen::Vector3d> const& points,
    Eigen::Vector3d const& origin,
    Eigen::Vector3d const& scaling_coefficients);

} // namespace geometry
} // namespace sbs

#endif // SBS_GEOMETRY_TRANSFORM_H
