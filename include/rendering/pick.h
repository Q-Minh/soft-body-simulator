#ifndef SBS_RENDERING_PICK_H
#define SBS_RENDERING_PICK_H

#include "common/mesh.h"
#include "common/primitive.h"

namespace sbs {
namespace rendering {

sbs::common::point_t unproject(
    Eigen::Vector3d const& coords,
    Eigen::Vector4d const& viewport,
    Eigen::Matrix4d const& projection,
    Eigen::Matrix4d const& view);

sbs::common::ray_t unproject_ray(
    Eigen::Vector2d const& coords,
    Eigen::Vector4d const& viewport,
    Eigen::Matrix4d const& projection,
    Eigen::Matrix4d const& view);

std::optional<
    std::tuple<std::uint32_t /* hit triangle */, double /* u */, double /* v */, double /* w */>>
pick(common::ray_t const& ray, common::shared_vertex_surface_mesh_t const& mesh);

} // namespace rendering
} // namespace sbs

#endif // SBS_RENDERING_PICK_H