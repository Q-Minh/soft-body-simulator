#ifndef SBS_RENDERING_PICK_H
#define SBS_RENDERING_PICK_H

#include "sbs/common/primitive.h"

#include <functional>
#include <memory>
#include <vector>

// forward declares
namespace sbs {
namespace common {

class shared_vertex_surface_mesh_i;

} // namespace common
} // namespace sbs

struct GLFWwindow;

namespace sbs {
namespace rendering {

// forward declare
class renderer_t;

class picker_t
{
  public:
    picker_t(
        renderer_t const* renderer,
        std::vector<std::shared_ptr<common::shared_vertex_surface_mesh_i>> const& nodes);

    void mouse_button_pressed_event(GLFWwindow* window, int button, int action, int mods);
    void mouse_moved_event(GLFWwindow* window, double x, double y);
    bool is_picking() const;
    bool is_usable() const;

    /**
     * Mandatory callbacks
     */
    std::function<bool(int /*button*/, int /*action*/, int /*mods*/)> should_picking_start;
    std::function<bool(
        bool /*left mouse button pressed*/,
        bool /*right mouse button pressed*/,
        bool /*middle mouse button pressed*/,
        bool /*alt key pressed*/,
        bool /*ctrl key pressed*/,
        bool /*shift key pressed*/)>
        should_picking_stop;
    std::function<bool(std::shared_ptr<common::shared_vertex_surface_mesh_i>)> should_pick;

    /**
     * Optional callbacks
     */
    std::function<void(std::shared_ptr<common::shared_vertex_surface_mesh_i>, std::uint32_t)>
        picked;
    std::function<void(std::shared_ptr<common::shared_vertex_surface_mesh_i>, std::uint32_t)>
        unpicked;
    std::function<void(
        Eigen::Vector3d const& /*d*/,
        std::shared_ptr<common::shared_vertex_surface_mesh_i> /*node*/,
        std::uint32_t /*picked vertex*/)>
        on_picker_moved;

  private:
    void unpick();

    renderer_t const* renderer_;
    std::vector<std::shared_ptr<common::shared_vertex_surface_mesh_i>> nodes_;
    std::vector<bool> is_node_picked_;
    std::vector<std::uint32_t> picked_vertices_;
    bool is_picking_;
    double xprev_, yprev_;
};

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