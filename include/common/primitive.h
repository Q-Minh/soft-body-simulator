#ifndef SBS_COMMON_PRIMITIVE_H
#define SBS_COMMON_PRIMITIVE_H

#include <Eigen/Core>
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
    point_t p, q;
};

struct triangle_t
{
    triangle_t(point_t const& a, point_t const& b, point_t const& c);

    normal_t normal() const;

    point_t a, b, c;
};

struct ray_t
{
    ray_t(point_t const& p, direction_t const& v, double t);

    point_t p;
    direction_t v;
    double t;
};

struct sphere_t;
struct aabb_t;

struct bounding_volume_t
{
    virtual bool contains(point_t const& point) const                              = 0;
    virtual bool contains(triangle_t const& triangle) const                        = 0;
    virtual std::optional<point_t> intersects(ray_t const& ray) const              = 0;
    virtual std::optional<point_t> intersects(line_segment_t const& segment) const = 0;
    virtual std::optional<point_t> intersects(triangle_t const& triangle) const    = 0;
    virtual std::optional<point_t> intersects(sphere_t const& sphere) const        = 0;
    virtual std::optional<point_t> intersects(aabb_t const& aabb) const            = 0;
};

struct sphere_t : public bounding_volume_t
{
    sphere_t(point_t const& center, double radius);

    point_t center;
    double radius;
};

struct aabb_t : public bounding_volume_t
{
    aabb_t(point_t const& center, vector3d_t const& extent);

    point_t center;
    vector3d_t extent;
};

std::tuple<double, double, double>
barycentric_coordinates(point_t const& A, point_t const& B, point_t const& C, point_t const& p);

std::optional<point_t> intersect(line_segment_t const& segment, triangle_t const& triangle);

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_PRIMITIVE_H