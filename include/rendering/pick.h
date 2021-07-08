#ifndef SBS_RENDERING_PICK_H
#define SBS_RENDERING_PICK_H

#include "common/primitive.h"

// forward declares
namespace sbs {
namespace common {

class shared_vertex_surface_mesh_i;

} // namespace common
} // namespace sbs

namespace sbs {
namespace rendering {

common::point_t unproject(
    Eigen::Vector3d const& coords,
    Eigen::Vector4d const& viewport,
    Eigen::Matrix4d const& projection,
    Eigen::Matrix4d const& view);

common::ray_t unproject_ray(
    Eigen::Vector2d const& coords,
    Eigen::Vector4d const& viewport,
    Eigen::Matrix4d const& projection,
    Eigen::Matrix4d const& view);

std::optional<
    std::tuple<std::uint32_t /* hit triangle */, double /* u */, double /* v */, double /* w */>>
pick(common::ray_t const& ray, common::shared_vertex_surface_mesh_i const& mesh);

std::optional<std::uint32_t>
pick_vertex(common::ray_t const& ray, common::shared_vertex_surface_mesh_i const& mesh);

} // namespace rendering
} // namespace sbs

#endif // SBS_RENDERING_PICK_H