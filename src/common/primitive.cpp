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

std::array<line_segment_t, 3u> triangle_t::edges() const
{
    return std::array<line_segment_t, 3u>{
        line_segment_t{p1(), p2()},
        line_segment_t{p2(), p3()},
        line_segment_t{p3(), p1()}};
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

sphere_t sphere_t::from(tetrahedron_t const& t)
{
    common::point_t const approx_barycenter =
        0.25 * t.p1() + 0.25 * t.p2() + 0.25 * t.p3() + 0.25 * t.p4();
    auto const d1     = (approx_barycenter - t.p1()).norm();
    auto const d2     = (approx_barycenter - t.p2()).norm();
    auto const d3     = (approx_barycenter - t.p3()).norm();
    auto const d4     = (approx_barycenter - t.p4()).norm();
    auto const radius = std::max({d1, d2, d3, d4});
    sphere_t const sphere{approx_barycenter, radius};
    return sphere;
}

sphere_t sphere_t::from(triangle_t const& t)
{
    common::point_t const approx_barycenter = 0.33 * t.p1() + 0.33 * t.p2() + 0.34 * t.p3();
    auto const d1                           = (approx_barycenter - t.p1()).norm();
    auto const d2                           = (approx_barycenter - t.p2()).norm();
    auto const d3                           = (approx_barycenter - t.p3()).norm();
    auto const radius                       = std::max({d1, d2, d3});
    sphere_t const sphere{approx_barycenter, radius};
    return sphere;
}

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

bool intersects(tetrahedron_t const& t1, tetrahedron_t const& t2)
{
    std::array<common::normal_t, 44u> separating_axis{};
    auto const& t1_triangles = t1.faces();
    auto const& t2_triangles = t2.faces();

    auto const face_axis_transform = [](common::triangle_t const& f) {
        return f.normal();
    };

    auto separating_axis_it = std::transform(
        t1_triangles.begin(),
        t1_triangles.end(),
        separating_axis.begin(),
        face_axis_transform);
    separating_axis_it = std::transform(
        t2_triangles.begin(),
        t2_triangles.end(),
        separating_axis_it,
        face_axis_transform);

    auto const t1_edges = t1.edges();
    auto const t2_edges = t2.edges();

    for (auto const& edge : t1_edges)
    {
        auto const edge_edge_separating_axis_transform_op =
            [&edge](common::line_segment_t const& e) {
                Eigen::Vector3d const d1 = edge.q - edge.p;
                Eigen::Vector3d d2       = e.q - e.p;

                auto axis = d1.cross(d2);

                double constexpr eps = std::numeric_limits<double>::epsilon();
                // Check if edges are parallel up to numerical precision eps
                if (axis.isZero(eps))
                {
                    d2   = e.p - edge.p;
                    axis = d1.cross(d2);
                }
                // Check if edges are on the same line up to numerical precision eps
                if (axis.isZero(eps))
                {
                    // Cancel this axis as a potential separating axis. When projecting
                    // nodes onto this axis, everything will be zeroed out, and thus no
                    // separating interval can be found. The separating axis test for
                    // this axis will fail.
                    axis.setZero();
                }

                return axis;
            };
        separating_axis_it = std::transform(
            t2_edges.begin(),
            t2_edges.end(),
            separating_axis_it,
            edge_edge_separating_axis_transform_op);
    }

    auto const is_separating_axis = [&t1, &t2](common::normal_t const& axis) {
        if (axis.isZero())
            return false;

        auto const project = [axis](common::point_t const& p) {
            return p.dot(axis);
        };

        std::array<double, 4u> const projection1{
            project(t1.p1()),
            project(t1.p2()),
            project(t1.p3()),
            project(t1.p4())};

        std::array<double, 4u> const projection2{
            project(t2.p1()),
            project(t2.p2()),
            project(t2.p3()),
            project(t2.p4())};

        auto const [min_it1, max_it1] = std::minmax_element(projection1.begin(), projection1.end());
        auto const [min_it2, max_it2] = std::minmax_element(projection2.begin(), projection2.end());

        bool const has_separating_interval = (*max_it1 < *min_it2) || (*max_it2 < *min_it1);
        return has_separating_interval;
    };

    return std::none_of(separating_axis.begin(), separating_axis.end(), is_separating_axis);
}

bool intersects(triangle_t const& triangle, tetrahedron_t const& tetrahedron)
{
    std::array<common::normal_t, 23u> separating_axis{};
    auto const& tet_triangles = tetrahedron.faces();

    auto const face_axis_transform = [](common::triangle_t const& f) {
        return f.normal();
    };

    separating_axis.front() = triangle.normal();

    auto separating_axis_it = separating_axis.begin() + 1u;
    separating_axis_it      = std::transform(
        tet_triangles.begin(),
        tet_triangles.end(),
        separating_axis_it,
        face_axis_transform);

    auto const t1_edges = triangle.edges();
    auto const t2_edges = tetrahedron.edges();

    for (auto const& edge : t1_edges)
    {
        auto const edge_edge_separating_axis_transform_op =
            [&edge](common::line_segment_t const& e) {
                Eigen::Vector3d const d1 = edge.q - edge.p;
                Eigen::Vector3d d2       = e.q - e.p;

                auto axis = d1.cross(d2);

                double constexpr eps = std::numeric_limits<double>::epsilon();
                // Check if edges are parallel up to numerical precision eps
                if (axis.isZero(eps))
                {
                    d2   = e.p - edge.p;
                    axis = d1.cross(d2);
                }
                // Check if edges are on the same line up to numerical precision eps
                if (axis.isZero(eps))
                {
                    // Cancel this axis as a potential separating axis. When projecting
                    // nodes onto this axis, everything will be zeroed out, and thus no
                    // separating interval can be found. The separating axis test for
                    // this axis will fail.
                    axis.setZero();
                }

                return axis;
            };
        separating_axis_it = std::transform(
            t2_edges.begin(),
            t2_edges.end(),
            separating_axis_it,
            edge_edge_separating_axis_transform_op);
    }

    auto const is_separating_axis = [&triangle, &tetrahedron](common::normal_t const& axis) {
        if (axis.isZero())
            return false;

        auto const project = [axis](common::point_t const& p) {
            return p.dot(axis);
        };

        std::array<double, 3u> const projection1{
            project(triangle.p1()),
            project(triangle.p2()),
            project(triangle.p3())};

        std::array<double, 4u> const projection2{
            project(tetrahedron.p1()),
            project(tetrahedron.p2()),
            project(tetrahedron.p3()),
            project(tetrahedron.p4())};

        auto const [min_it1, max_it1] = std::minmax_element(projection1.begin(), projection1.end());
        auto const [min_it2, max_it2] = std::minmax_element(projection2.begin(), projection2.end());

        bool const has_separating_interval = (*max_it1 < *min_it2) || (*max_it2 < *min_it1);
        return has_separating_interval;
    };

    return std::none_of(separating_axis.begin(), separating_axis.end(), is_separating_axis);
}

bool intersects(point_t const& point, tetrahedron_t const& tetrahedron)
{
    auto const project = [](common::point_t const& p, common::triangle_t const& triangle) {
        auto const n = triangle.normal();
        auto const d = p - triangle.p1();
        return n.dot(d);
    };
    std::array<common::triangle_t, 4u> const faces = tetrahedron.faces();
    std::array<double, 4u> const projections{
        project(point, faces[0]),
        project(point, faces[1]),
        project(point, faces[2]),
        project(point, faces[3])};

    return std::none_of(projections.begin(), projections.end(), [](double const s) {
        return s > 0.;
    });
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

std::optional<point_t> intersect_twoway(line_segment_t const& segment, triangle_t const& triangle)
{
    auto const intersection = intersect(segment, triangle);
    if (intersection.has_value())
        return intersection;

    line_segment_t const flipped_segment{segment.q, segment.p};
    return intersect(flipped_segment, triangle);
}

/**
 * @brief
 * Implementation of closest point on triangle to a point P from Christer Ericson's Real-Time
 * Collision Detection
 * @param p Point off triangle
 * @param t Triangle on which we wish to find the closest point to p
 * @return The closest point q on triangle t to point p
 */
point_t closest_point(point_t const& p, triangle_t const& t)
{
    auto const& a = t.a();
    auto const& b = t.b();
    auto const& c = t.c();

    // Check if P in vertex region outside A
    common::vector3d_t const ab = b - a;
    common::vector3d_t const ac = c - a;
    common::vector3d_t const ap = p - a;
    double const d1             = ab.dot(ap);
    double const d2             = ac.dot(ap);
    if (d1 <= 0.0 && d2 <= 0.0)
        return a; // barycentric coordinates (1,0,0)

    // Check if P in vertex region outside B
    common::vector3d_t const bp = p - b;
    double const d3             = ab.dot(bp);
    double const d4             = ac.dot(bp);
    if (d3 >= 0.0 && d4 <= d3)
        return b; // barycentric coordinates (0,1,0)

    // Check if P in edge region of AB, if so return projection of P onto AB
    double const vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        double const v = d1 / (d1 - d3);
        return a + v * ab; // barycentric coordinates (1-v, v,0)
    }

    // Check if P in vertex region outside C
    common::vector3d_t const cp = p - c;
    double const d5             = ab.dot(cp);
    double const d6             = ac.dot(cp);
    if (d6 >= 0.0 && d5 <= d6)
        return c; // barycentric coordinates (0,0,1)

    // Check if P in edge region of AC, if so return projection of P onto AC
    double const vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        double const w = d2 / (d2 - d6);
        return a + w * ac; // barycentric coordinates (1-w, 0, w)
    }

    // Check if P in edge region of BC, if so return projection of P onto BC
    double const va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        double const w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + w * (c - b); // barycentric coordinates (0,1-w, w)
    }

    // P inside face region. Compute Q through its barycentric coordinates (u, v, w)
    double const denom = 1.0 / (va + vb + vc);
    double const v     = vb * denom;
    double const w     = vc * denom;
    return a + ab * v + ac * w; //=u*a+v*b+w*c,u=va* denom=1.0f−v−w
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