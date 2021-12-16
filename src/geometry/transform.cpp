#include "sbs/geometry/transform.h"

#include "..\..\include\sbs\geometry\transform.h"

#include <Eigen/Geometry>

namespace sbs {
namespace geometry {
std::vector<Eigen::Vector3d>
translate(std::vector<Eigen::Vector3d> const& points, Eigen::Vector3d const& translation)
{
    std::vector<Eigen::Vector3d> translated_points{};
    translated_points.reserve(points.size());
    for (auto const& p : points)
    {
        Eigen::Vector3d const translated_point = p + translation;
        translated_points.push_back(translated_point);
    }
    return translated_points;
}
Eigen::Vector3d center_of_geometry(std::vector<Eigen::Vector3d> const& points)
{
    scalar_type const N = static_cast<scalar_type>(points.size());
    Eigen::Vector3d sum{0., 0., 0.};
    for (auto const& p : points)
    {
        sum += p;
    }
    Eigen::Vector3d const mean = sum / N;
    return mean;
}

std::vector<Eigen::Vector3d> rotate(
    std::vector<Eigen::Vector3d> const& points,
    Eigen::AngleAxis<scalar_type> const& angle_axis,
    Eigen::Vector3d const& origin)
{
    Eigen::Matrix3d const rotation = angle_axis.toRotationMatrix();

    std::vector<Eigen::Vector3d> rotated_points{};
    rotated_points.reserve(points.size());

    for (Eigen::Vector3d const& p : points)
    {
        Eigen::Vector3d const centered_point = p - origin;
        // rotate around center of geometry, then translate back to initial frame
        Eigen::Vector3d const rotated_point = rotation * (centered_point) + origin;
        rotated_points.push_back(rotated_point);
    }

    return rotated_points;
}

std::vector<Eigen::Vector3d> scale(
    std::vector<Eigen::Vector3d> const& points,
    Eigen::Vector3d const& origin,
    Eigen::Vector3d const& scaling_coefficients)
{
    std::vector<Eigen::Vector3d> scaled_points{};
    scaled_points.reserve(points.size());

    auto const& sx = scaling_coefficients.x();
    auto const& sy = scaling_coefficients.y();
    auto const& sz = scaling_coefficients.z();

    for (Eigen::Vector3d const& p : points)
    {
        Eigen::Vector3d const scaled_point{sx * p.x(), sy * p.y(), sz * p.z()};
        scaled_points.push_back(scaled_point);
    }
    return scaled_points;
}

} // namespace geometry
} // namespace sbs