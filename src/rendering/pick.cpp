#include "sbs/rendering/pick.h"

#include "sbs/common/mesh.h"
#include "sbs/rendering/renderer.h"

#include <Eigen/LU>
#include <GLFW/glfw3.h>
#include <iostream>

namespace sbs {
namespace rendering {

picker_t::picker_t(
    renderer_t const* renderer,
    std::vector<std::shared_ptr<common::shared_vertex_surface_mesh_i>> const& nodes)
    : renderer_(renderer),
      nodes_(nodes),
      is_node_picked_(nodes.size(), false),
      picked_vertices_(nodes.size(), 0u),
      is_picking_{false},
      xprev_{0.},
      yprev_{0.}
{
}

void picker_t::mouse_button_pressed_event(GLFWwindow* window, int button, int action, int mods)
{
    if (!should_picking_start)
        return;

    is_picking_ = should_picking_start(button, action, mods);
    if (is_picking_)
    {
        if (!should_pick)
            return;

        int viewport_gl[4];
        glGetIntegerv(GL_VIEWPORT, viewport_gl);
        Eigen::Vector4d const viewport{
            static_cast<double>(viewport_gl[0]),
            static_cast<double>(viewport_gl[1]),
            static_cast<double>(viewport_gl[2]),
            static_cast<double>(viewport_gl[3])};

        float const aspect_ratio =
            static_cast<float>(viewport(2)) / static_cast<float>(viewport(3));
        Eigen::Matrix4d const projection = renderer_->camera().projection_matrix(aspect_ratio);
        Eigen::Matrix4d const view       = renderer_->camera().view_matrix();

        double x, y;
        glfwGetCursorPos(window, &x, &y);
        xprev_ = x;
        yprev_ = y;

        auto const ray = sbs::rendering::unproject_ray({x, y}, viewport, projection, view);

        for (std::size_t i = 0u; i < nodes_.size(); ++i)
        {
            auto const& node = nodes_[i];

            if (!should_pick(node))
                continue;

            auto const vid = rendering::pick_vertex(ray, *node);
            if (!vid.has_value())
                continue;

            is_node_picked_[i]  = true;
            picked_vertices_[i] = *vid;

            if (!picked)
                continue;

            picked(node, *vid);
        }
    }
    else
    {
        unpick();
    }
}

void picker_t::mouse_moved_event(GLFWwindow* window, double x, double y)
{
    if (!should_picking_stop)
    {
        return;
    }

    if (!is_picking_)
    {
        return;
    }

    bool const left_mouse_button_pressed =
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    bool const right_mouse_button_pressed =
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
    bool const middle_mouse_button_pressed =
        glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;

    bool const alt_key_pressed = glfwGetKey(window, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
                                 glfwGetKey(window, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
    bool const ctrl_key_pressed = glfwGetKey(window, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
                                  glfwGetKey(window, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;
    bool const shift_key_pressed = glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
                                   glfwGetKey(window, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;

    is_picking_ = !should_picking_stop(
        left_mouse_button_pressed,
        right_mouse_button_pressed,
        middle_mouse_button_pressed,
        alt_key_pressed,
        ctrl_key_pressed,
        shift_key_pressed);

    if (!is_picking_)
    {
        unpick();
    }
    else
    {
        int viewport_gl[4];
        glGetIntegerv(GL_VIEWPORT, viewport_gl);
        Eigen::Vector4d const viewport{
            static_cast<double>(viewport_gl[0]),
            static_cast<double>(viewport_gl[1]),
            static_cast<double>(viewport_gl[2]),
            static_cast<double>(viewport_gl[3])};

        float const aspect_ratio =
            static_cast<float>(viewport(2)) / static_cast<float>(viewport(3));
        Eigen::Matrix4d const projection = renderer_->camera().projection_matrix(aspect_ratio);
        Eigen::Matrix4d const view       = renderer_->camera().view_matrix();

        for (std::size_t i = 0u; i < nodes_.size(); ++i)
        {
            if (!is_node_picked_[i])
                continue;

            if (!on_picker_moved)
                continue;

            common::point_t const p1 = rendering::unproject(
                Eigen::Vector3d{xprev_, yprev_, 0.},
                viewport,
                projection,
                view);
            common::point_t const p2 =
                rendering::unproject(Eigen::Vector3d{x, y, 0.}, viewport, projection, view);
            Eigen::Vector3d const d{p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z()};

            on_picker_moved(d, nodes_[i], picked_vertices_[i]);
        }
    }

    xprev_ = x;
    yprev_ = y;
}

bool picker_t::is_picking() const
{
    return is_picking_;
}

bool picker_t::is_usable() const
{
    return static_cast<bool>(should_picking_start) && static_cast<bool>(should_picking_stop) &&
           static_cast<bool>(should_pick);
}

void picker_t::unpick()
{
    for (std::size_t i = 0u; i < nodes_.size(); ++i)
    {
        if (is_node_picked_[i])
        {
            is_node_picked_[i] = false;
            if (unpicked)
            {
                unpicked(nodes_[i], picked_vertices_[i]);
            }
        }
    }
}

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
pick(common::ray_t const& ray, common::shared_vertex_surface_mesh_i const& mesh)
{
    using return_type = std::optional<std::tuple<std::uint32_t, double, double, double>>;

    std::size_t const num_triangles = mesh.triangle_count();
    for (std::size_t f = 0u; f < num_triangles; ++f)
    {
        auto const t  = mesh.triangle(f);
        auto const v1 = mesh.vertex(t.v1);
        auto const v2 = mesh.vertex(t.v2);
        auto const v3 = mesh.vertex(t.v3);

        auto const& a = common::point_t{v1.x, v1.y, v1.z};
        auto const& b = common::point_t{v2.x, v2.y, v2.z};
        auto const& c = common::point_t{v3.x, v3.y, v3.z};

        common::triangle_t const triangle{a, b, c};
        auto const intersection = common::intersect(ray, triangle);

        if (!intersection.has_value())
            continue;

        auto const [u, v, w] = common::barycentric_coordinates(a, b, c, intersection.value());
        return std::make_tuple(static_cast<std::uint32_t>(f), u, v, w);
    }

    return {};
}
std::optional<std::uint32_t>
pick_vertex(common::ray_t const& ray, common::shared_vertex_surface_mesh_i const& mesh)
{
    auto const picked_triangle = pick(ray, mesh);
    if (!picked_triangle.has_value())
        return {};

    auto const& [f, u, v, w] = picked_triangle.value();
    std::uint32_t vi         = mesh.triangle(f).v1;
    if (v > u && v > w)
        vi = mesh.triangle(f).v2;
    if (w > u && w > v)
        vi = mesh.triangle(f).v3;

    return vi;
};

} // namespace rendering
} // namespace sbs