#include "common/primitive.h"

namespace sbs {
namespace common {

line_segment_t::line_segment_t(point_t const& p, point_t const& q) : p(p), q(q) {}

triangle_t::triangle_t(point_t const& a, point_t const& b, point_t const& c) : a(a), b(b), c(c) {}

normal_t triangle_t::normal() const
{
    Eigen::Vector3d const ab = b - a;
    Eigen::Vector3d const ac = c - a;
    return ab.cross(ac).normalized();
}

ray_t::ray_t(point_t const& p, direction_t const& v, double t) : p(p), v(v), t(t) {}

sphere_t::sphere_t(point_t const& center, double radius) : center(center), radius(radius) {}

aabb_t::aabb_t(point_t const& center, vector3d_t const& extent) : center(center), extent(extent) {}

std::tuple<double, double, double>
barycentric_coordinates(point_t const& A, point_t const& B, point_t const& C, point_t const& p)
{
    Eigen::Matrix3d M{};
    M.col(0u) = A;
    M.col(1u) = B;
    M.col(2u) = C;

    Eigen::Vector3d const uvw = M.inverse() * p;
    return std::make_tuple(uvw(0u), uvw(1u), uvw(2u));
}

} // namespace common
} // namespace sbs