#ifndef SBS_COMMON_PRIMITIVE_H
#define SBS_COMMON_PRIMITIVE_H

#include <Eigen/Core>
#include <array>
#include <optional>

namespace sbs {
namespace common {

using point_t     = Eigen::Vector3d;
using normal_t    = Eigen::Vector3d;
using direction_t = Eigen::Vector3d;
using vector3d_t  = Eigen::Vector3d;

struct line_segment_t
{
    line_segment_t(point_t const& p, point_t const& q);

    friend bool operator==(line_segment_t const& l1, line_segment_t const& l2);
    friend bool operator!=(line_segment_t const& l1, line_segment_t const& l2);
    friend line_segment_t operator*(Eigen::Matrix3d const& R, line_segment_t const& l);
    friend line_segment_t operator+(line_segment_t const& l, vector3d_t const& t);
    friend line_segment_t operator+(vector3d_t const& t, line_segment_t const& l);

    point_t p, q;
};

struct triangle_t
{
    triangle_t(point_t const& a, point_t const& b, point_t const& c);

    normal_t normal() const;
    double area() const;

    std::array<line_segment_t, 3u> edges() const;

    std::array<point_t, 3u> const& nodes() const;
    std::array<point_t, 3u>& nodes();

    point_t const& p1() const;
    point_t const& p2() const;
    point_t const& p3() const;

    point_t& p1();
    point_t& p2();
    point_t& p3();

    point_t const& a() const;
    point_t const& b() const;
    point_t const& c() const;

    point_t& a();
    point_t& b();
    point_t& c();

  private:
    std::array<point_t, 3u> p_;
};

struct tetrahedron_t
{
    tetrahedron_t(point_t const& p1, point_t const& p2, point_t const& p3, point_t const& p4);

    double unsigned_volume() const;
    double signed_volume() const;
    std::array<triangle_t, 4u> faces() const;
    std::array<line_segment_t, 6u> edges() const;
    std::array<point_t, 4u> const& nodes() const;
    std::array<point_t, 4u>& nodes();

    point_t const& p1() const;
    point_t const& p2() const;
    point_t const& p3() const;
    point_t const& p4() const;

    point_t& p1();
    point_t& p2();
    point_t& p3();
    point_t& p4();

  private:
    std::array<point_t, 4u> p_;
};

struct ray_t
{
    ray_t(point_t const& p, direction_t const& v);

    point_t p;
    direction_t v;
};

struct plane_t
{
    plane_t(point_t const& p, normal_t const& n);
    plane_t(triangle_t const& t);

    double signed_distance(point_t const& q) const;

    point_t p;
    normal_t n;
};

struct sphere_t;
struct aabb_t;

struct bounding_volume_t
{
    // virtual bool contains(point_t const& point) const                              = 0;
    // virtual bool contains(triangle_t const& triangle) const                        = 0;
    // virtual std::optional<point_t> intersects(ray_t const& ray) const              = 0;
    // virtual std::optional<point_t> intersects(line_segment_t const& segment) const = 0;
    // virtual std::optional<point_t> intersects(triangle_t const& triangle) const    = 0;
    // virtual std::optional<point_t> intersects(sphere_t const& sphere) const        = 0;
    // virtual std::optional<point_t> intersects(aabb_t const& aabb) const            = 0;
};

struct sphere_t : public bounding_volume_t
{
    sphere_t(point_t const& center, double radius);

    static sphere_t from(tetrahedron_t const& t);
    static sphere_t from(triangle_t const& t);

    point_t center;
    double radius;
};

struct aabb_t : public bounding_volume_t
{
    aabb_t(point_t const& min, point_t const& max);

    bool contains(point_t const& p) const;
    static aabb_t from(tetrahedron_t const& t);
    static aabb_t from(triangle_t const& t);

    point_t min, max;
};

std::tuple<double, double, double>
barycentric_coordinates(point_t const& A, point_t const& B, point_t const& C, point_t const& p);

bool intersects(tetrahedron_t const& t1, tetrahedron_t const& t2);
bool intersects(triangle_t const& triangle, tetrahedron_t const& tetrahedron);
bool intersects(point_t const& p, tetrahedron_t const& tetrahedron);

std::optional<point_t> intersect(line_segment_t const& segment, triangle_t const& triangle);
std::optional<point_t> intersect(ray_t const& ray, triangle_t const& triangle);
std::optional<point_t> intersect_twoway(line_segment_t const& segment, triangle_t const& triangle);
std::optional<point_t> intersect_twoway(ray_t const& ray, triangle_t const& triangle);

point_t closest_point(point_t const& p, triangle_t const& t);

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_PRIMITIVE_H