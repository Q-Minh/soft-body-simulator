#include "rendering/pick.h"

#include <Eigen/LU>

namespace sbs {
namespace rendering {

sbs::common::point_t unproject(
    Eigen::Vector3d const& coords,
    Eigen::Vector4d const& viewport,
    Eigen::Matrix4d const& projection,
    Eigen::Matrix4d const& view)
{
    double const x = (coords.x() - viewport(0)) / viewport(2);
    double const y = 1. - (coords.y() - viewport(1)) / viewport(3);

    Eigen::Vector4d const projective_screen_space_point =
        Eigen::Vector4d{x, y, coords.z(), 1.}.array() * 2. - 1.;

    Eigen::Matrix4d const screen_to_world_transform = (projection * view).inverse();

    Eigen::Vector4d const projective_world_space_point =
        screen_to_world_transform * projective_screen_space_point;

    Eigen::Vector3d const world_space_point =
        (projective_world_space_point / projective_world_space_point.w()).head(3);

    return world_space_point;
};

sbs::common::ray_t unproject_ray(
    Eigen::Vector2d const& coords,
    Eigen::Vector4d const& viewport,
    Eigen::Matrix4d const& projection,
    Eigen::Matrix4d const& view)
{
    Eigen::Vector3d const win_source(coords.x(), coords.y(), 0.);
    Eigen::Vector3d const win_dest(coords.x(), coords.y(), 1.);

    Eigen::Vector3d const source = unproject(win_source, viewport, projection, view);
    Eigen::Vector3d const dest   = unproject(win_dest, viewport, projection, view);
    Eigen::Vector3d const dir    = dest - source;

    sbs::common::ray_t const ray{source, dir};
    return ray;
}

std::optional<std::tuple<std::uint32_t, double, double, double>>
pick(common::ray_t const& ray, common::shared_vertex_surface_mesh_t const& mesh)
{
    using return_type = std::optional<std::tuple<std::uint32_t, double, double, double>>;

    std::size_t const num_triangles = static_cast<std::size_t>(mesh.triangles().cols());
    for (std::size_t f = 0u; f < num_triangles; ++f)
    {
        auto const v1 = mesh.triangles()(0u, f);
        auto const v2 = mesh.triangles()(1u, f);
        auto const v3 = mesh.triangles()(2u, f);

        auto const& a = mesh.vertices().col(v1);
        auto const& b = mesh.vertices().col(v2);
        auto const& c = mesh.vertices().col(v3);

        common::triangle_t const triangle{a, b, c};
        auto const intersection = common::intersect(ray, triangle);

        if (!intersection.has_value())
            continue;

        auto const [u, v, w] = common::barycentric_coordinates(a, b, c, intersection.value());
        return std::make_tuple(static_cast<std::uint32_t>(f), u, v, w);
    }

    return {};
};

} // namespace rendering
} // namespace sbs