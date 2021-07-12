#include "common/primitive.h"

#include <Eigen/Geometry>

namespace sbs {
namespace common {

line_segment_t::line_segment_t(point_t const& p, point_t const& q) : p(p), q(q) {}

triangle_t::triangle_t(point_t const& a, point_t const& b, point_t const& c) : p_{a, b, c} {}

normal_t triangle_t::normal() const
{
    Eigen::Vector3d const ab = b() - a();
    Eigen::Vector3d const ac = c() - a();
    return ab.cross(ac).normalized();
}

double triangle_t::area() const
{
    return 0.5 * (b() - a()).cross(c() - a()).norm();
}

std::array<point_t, 3u> const& triangle_t::nodes() const
{
    return p_;
}

std::array<point_t, 3u>& triangle_t::nodes()
{
    return p_;
}

point_t const& triangle_t::p1() const
{
    return p_[0];
}

point_t const& triangle_t::p2() const
{
    return p_[1];
}

point_t const& triangle_t::p3() const
{
    return p_[2];
}

point_t& triangle_t::p1()
{
    return p_[0];
}

point_t& triangle_t::p2()
{
    return p_[1];
}

point_t& triangle_t::p3()
{
    return p_[2];
}

point_t const& triangle_t::a() const
{
    return p1();
}

point_t const& triangle_t::b() const
{
    return p2();
}

point_t const& triangle_t::c() const
{
    return p3();
}

point_t& triangle_t::a()
{
    return p1();
}

point_t& triangle_t::b()
{
    return p2();
}

point_t& triangle_t::c()
{
    return p3();
}

ray_t::ray_t(point_t const& p, direction_t const& v) : p(p), v(v) {}

sphere_t::sphere_t(point_t const& center, double radius) : center(center), radius(radius) {}

bool operator==(line_segment_t const& l1, line_segment_t const& l2)
{
    return l1.p.isApprox(l2.p) && l1.q.isApprox(l2.q);
}

bool operator!=(line_segment_t const& l1, line_segment_t const& l2)
{
    return !(l1 == l2);
}

line_segment_t operator+(line_segment_t const& l, vector3d_t const& t)
{
    return line_segment_t(l.p + t, l.q + t);
}

line_segment_t operator+(vector3d_t const& t, line_segment_t const& l)
{
    return l + t;
}

line_segment_t operator*(Eigen::Matrix3d const& R, line_segment_t const& l)
{
    return line_segment_t(R * l.p, R * l.q);
}

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
    Eigen::Vector3d const ab = triangle.b() - triangle.a();
    Eigen::Vector3d const ac = triangle.c() - triangle.a();
    Eigen::Vector3d const qp = segment.p - segment.q;

    Eigen::Vector3d const n = ab.cross(ac);

    double const d = qp.dot(n);
    if (d <= 0.)
        return {};

    Eigen::Vector3d const ap = segment.p - triangle.a();
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
    point_t const intersection = u * triangle.a() + v * triangle.b() + w * triangle.c();
    return intersection;
}

std::optional<point_t> intersect(ray_t const& ray, triangle_t const& triangle)
{
    Eigen::Vector3d const& p = ray.p;
    Eigen::Vector3d const& q = ray.p + 1. * ray.v;
    Eigen::Vector3d const ab = triangle.b() - triangle.a();
    Eigen::Vector3d const ac = triangle.c() - triangle.a();
    Eigen::Vector3d const qp = p - q;

    Eigen::Vector3d const n = ab.cross(ac);

    double const d = qp.dot(n);
    if (d <= 0.)
        return {};

    Eigen::Vector3d const ap = p - triangle.a();
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
    point_t const intersection = u * triangle.a() + v * triangle.b() + w * triangle.c();
    return intersection;
}

tetrahedron_t::tetrahedron_t(
    point_t const& p1,
    point_t const& p2,
    point_t const& p3,
    point_t const& p4)
    : p_{p1, p2, p3, p4}
{
}

double tetrahedron_t::unsigned_volume() const
{
    return std::abs(signed_volume());
}

double tetrahedron_t::signed_volume() const
{
    vector3d_t const p21 = p2() - p1();
    vector3d_t const p31 = p3() - p1();
    vector3d_t const p41 = p4() - p1();
    return p21.cross(p31).dot(p41);
}

std::array<triangle_t, 4u> tetrahedron_t::faces() const
{
    return std::array<triangle_t, 4u>{
        triangle_t{p1(), p2(), p4()},
        triangle_t{p2(), p3(), p4()},
        triangle_t{p3(), p1(), p4()},
        triangle_t{p1(), p3(), p2()}};
}

std::array<line_segment_t, 6u> tetrahedron_t::edges() const
{
    return std::array<line_segment_t, 6u>{
        line_segment_t{p1(), p2()},
        line_segment_t{p2(), p3()},
        line_segment_t{p3(), p1()},
        line_segment_t{p1(), p4()},
        line_segment_t{p2(), p4()},
        line_segment_t{p3(), p4()},
    };
}

std::array<point_t, 4u> const& tetrahedron_t::nodes() const
{
    return p_;
}

std::array<point_t, 4u>& tetrahedron_t::nodes()
{
    return p_;
}

point_t const& tetrahedron_t::p1() const
{
    return p_[0];
}

point_t const& tetrahedron_t::p2() const
{
    return p_[1];
}

point_t const& tetrahedron_t::p3() const
{
    return p_[2];
}

point_t const& tetrahedron_t::p4() const
{
    return p_[3];
}

point_t& tetrahedron_t::p1()
{
    return p_[0];
}

point_t& tetrahedron_t::p2()
{
    return p_[1];
}

point_t& tetrahedron_t::p3()
{
    return p_[2];
}

point_t& tetrahedron_t::p4()
{
    return p_[3];
}

aabb_t::aabb_t(point_t const& min, point_t const& max) : min(min), max(max) {}

bool aabb_t::contains(point_t const& p) const
{
    return (p.x() >= min.x() && p.x() <= max.x()) && (p.y() >= min.y() && p.y() <= max.y()) &&
           (p.z() >= min.z() && p.z() <= max.z());
}

aabb_t aabb_t::from(tetrahedron_t const& t)
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    point_t min{inf, inf, inf}, max{-inf, -inf, -inf};
    for (auto const& p : t.nodes())
    {
        if (p.x() > max.x())
            max.x() = p.x();
        if (p.y() > max.y())
            max.y() = p.y();
        if (p.z() > max.z())
            max.z() = p.z();

        if (p.x() < min.x())
            min.x() = p.x();
        if (p.y() < min.y())
            min.y() = p.y();
        if (p.z() < min.z())
            min.z() = p.z();
    }
    return aabb_t{min, max};
}

aabb_t aabb_t::from(triangle_t const& t)
{
    constexpr double inf = std::numeric_limits<double>::infinity();
    point_t min{inf, inf, inf}, max{-inf, -inf, -inf};
    for (auto const& p : t.nodes())
    {
        if (p.x() > max.x())
            max.x() = p.x();
        if (p.y() > max.y())
            max.y() = p.y();
        if (p.z() > max.z())
            max.z() = p.z();

        if (p.x() < min.x())
            min.x() = p.x();
        if (p.y() < min.y())
            min.y() = p.y();
        if (p.z() < min.z())
            min.z() = p.z();
    }
    return aabb_t{min, max};
}

} // namespace common
} // namespace sbs