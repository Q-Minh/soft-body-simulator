#ifndef SBS_PHYSICS_COLLISION_INTERSECTIONS_H
#define SBS_PHYSICS_COLLISION_INTERSECTIONS_H

#include <Eigen/Core>
#include <optional>

namespace sbs {
namespace physics {
namespace collision {

using point_t     = Eigen::Vector3d;
using normal_t    = Eigen::Vector3d;
using direction_t = Eigen::Vector3d;

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

std::optional<point_t> intersect(line_segment_t const& segment, triangle_t const& triangle);

} // namespace collision
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_INTERSECTIONS_H