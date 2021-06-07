#include "common/primitive.h"

#include <Eigen/Geometry>

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

ray_t::ray_t(point_t const& p, direction_t const& v) : p(p), v(v) {}

sphere_t::sphere_t(point_t const& center, double radius) : center(center), radius(radius) {}

aabb_t::aabb_t(point_t const& center, vector3d_t const& extent) : center(center), extent(extent) {}

std::tuple<double, double, double>
barycentric_coordinates(point_t const& A, point_t const& B, point_t const& C, point_t const& p)
{
    Eigen::Vector3d const v0 = B - A;
    Eigen::Vector3d const v1 = C - A;
    Eigen::Vector3d const v2 = p - A;

    Eigen::Vector3d const AB = B - A;
    Eigen::Vector3d const AC = C - A;
    Eigen::Vector3d const AP = p - A;

    double const d00   = AB.dot(AB);
    double const d01   = AB.dot(AC);
    double const d11   = AC.dot(AC);
    double const d20   = AP.dot(AB);
    double const d21   = AP.dot(AC);
    double const denom = d00 * d11 - d01 * d01;
    double const v     = (d11 * d20 - d01 * d21) / denom;
    double const w     = (d00 * d21 - d01 * d20) / denom;
    double const u     = 1.0 - v - w;

    return std::make_tuple(u, v, w);
}

std::optional<point_t> intersect(line_segment_t const& segment, triangle_t const& triangle)
{
    Eigen::Vector3d const ab = triangle.b - triangle.a;
    Eigen::Vector3d const ac = triangle.c - triangle.a;
    Eigen::Vector3d const qp = segment.p - segment.q;

    Eigen::Vector3d const n = ab.cross(ac);

    double const d = qp.dot(n);
    if (d <= 0.)
        return {};

    Eigen::Vector3d const ap = segment.p - triangle.a;
    double const t           = ap.dot(n);
    if (t < 0.)
        return {};
    if (t > d)
        return {};

    Eigen::Vector3d const e = qp.cross(ap);
    double v                = ac.dot(e);
    if (v < 0. || v > d)
        return {};

    double w = -ab.dot(e);
    if (w < 0. || (v + w) > d)
        return {};

    double const ood = 1. / d;
    v *= ood;
    w *= ood;
    double const u             = 1. - v - w;
    point_t const intersection = u * triangle.a + v * triangle.b + w * triangle.c;
    return intersection;
}

std::optional<point_t> intersect(ray_t const& ray, triangle_t const& triangle)
{
    Eigen::Vector3d const& p = ray.p;
    Eigen::Vector3d const& q = ray.p + 1. * ray.v;
    Eigen::Vector3d const ab = triangle.b - triangle.a;
    Eigen::Vector3d const ac = triangle.c - triangle.a;
    Eigen::Vector3d const qp = p - q;

    Eigen::Vector3d const n = ab.cross(ac);

    double const d = qp.dot(n);
    if (d <= 0.)
        return {};

    Eigen::Vector3d const ap = p - triangle.a;
    double const t           = ap.dot(n);
    if (t < 0.)
        return {};

    Eigen::Vector3d const e = qp.cross(ap);
    double v                = ac.dot(e);
    if (v < 0. || v > d)
        return {};

    double w = -ab.dot(e);
    if (w < 0. || (v + w) > d)
        return {};

    double const ood = 1. / d;
    // t *= ood;
    v *= ood;
    w *= ood;
    double const u             = 1. - v - w;
    point_t const intersection = u * triangle.a + v * triangle.b + w * triangle.c;
    return intersection;
}

} // namespace common
} // namespace sbs