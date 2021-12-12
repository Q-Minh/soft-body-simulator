#ifndef SBS_RENDERING_RENDERER_H
#define SBS_RENDERING_RENDERER_H

// clang-format off
#include <glad/glad.h>
#include <GLFW/glfw3.h>
// clang-format on

#include "camera.h"
#include "light.h"
#include "sbs/aliases.h"
#include "sbs/common/node.h"
#include "shader.h"

#include <array>
#include <filesystem>
#include <functional>
#include <list>

namespace sbs {
namespace physics {
namespace xpbd {
// Forward declares
class simulation_t;

} // namespace xpbd
} // namespace physics

namespace rendering {

// forward declare
class picker_t;

class renderer_base_t
{
  public:
    virtual void framebuffer_size_callback(GLFWwindow* window, int width, int height)        = 0;
    virtual void mouse_move_callback(GLFWwindow* window, double xpos, double ypos)           = 0;
    virtual void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)         = 0;
    virtual void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) = 0;

    static renderer_base_t* active_renderer;

    virtual void set_as_active_renderer() { active_renderer = this; }

  protected:
    static void framebuffer_size_callback_dispatcher(GLFWwindow* window, int width, int height)
    {
        if (active_renderer != nullptr)
            active_renderer->framebuffer_size_callback(window, width, height);
    }
    static void mouse_move_callback_dispatcher(GLFWwindow* window, double xpos, double ypos)
    {
        if (active_renderer != nullptr)
            active_renderer->mouse_move_callback(window, xpos, ypos);
    }
    static void scroll_callback_dispatcher(GLFWwindow* window, double xoffset, double yoffset)
    {
        if (active_renderer != nullptr)
            active_renderer->scroll_callback(window, xoffset, yoffset);
    }
    static void
    mouse_button_callback_dispatcher(GLFWwindow* window, int button, int action, int mods)
    {
        if (active_renderer != nullptr)
            active_renderer->mouse_button_callback(window, button, action, mods);
    }
};

class renderer_t : public renderer_base_t
{
  public:
    renderer_t(
        physics::xpbd::simulation_t& simulation,
        point_light_t const& point_light,
        directional_light_t const& directional_light);

    ~renderer_t();

    bool initialize();

    bool use_mesh_shaders(
        std::filesystem::path const& vertex_shader_path,
        std::filesystem::path const& fragment_shader_path);

    bool use_wireframe_shaders(
        std::filesystem::path const& vertex_shader_path,
        std::filesystem::path const& fragment_shader_path);

    bool use_point_shaders(
        std::filesystem::path const& vertex_shader_path,
        std::filesystem::path const& fragment_shader_path);

    void launch();
    void close();

    std::vector<std::string> get_error_messages() const { return mesh_shader_.error_messages(); }

    std::uint32_t constexpr get_initial_window_width() const { return 800u; }
    std::uint32_t constexpr get_initial_window_height() const { return 600u; }

    camera_t const& camera() const { return camera_; }
    camera_t& camera() { return camera_; }

    void add_point(std::array<float, 9u> const& xyz_nxnynz_rgb_point);
    void clear_points();

    void add_rendered_object(std::unique_ptr<common::renderable_node_t> rendered_object);
    std::unique_ptr<common::renderable_node_t> const& rendered_object(index_type i) const
    {
        return rendered_objects_[i];
    }
    std::unique_ptr<common::renderable_node_t>& rendered_object(index_type i)
    {
        return rendered_objects_[i];
    }
    std::vector<std::unique_ptr<common::renderable_node_t>> const& rendered_objects() const
    {
        return rendered_objects_;
    }
    std::vector<std::unique_ptr<common::renderable_node_t>>& rendered_objects()
    {
        return rendered_objects_;
    }
    std::size_t rendered_object_count() const { return rendered_objects_.size(); }

    /**
     * Render loop hooks
     */
    std::function<void(physics::xpbd::simulation_t&)> on_new_imgui_frame;
    std::function<void(double /*render_frame_dt*/, physics::xpbd::simulation_t& /*scene*/)>
        on_new_physics_timestep;
    std::function<void(physics::xpbd::simulation_t&)> on_pre_render;

    /**
     * User input hooks
     */
    std::function<bool(GLFWwindow* /*window*/, double /*mouse_x_pos*/, double /*mouse_y_pos*/)>
        on_mouse_moved;
    std::function<void(GLFWwindow*, int, int, int)> on_mouse_button_pressed;

    /**
     * @brief Add picking features to the renderer
     */
    std::list<picker_t> pickers;

  protected:
    virtual void framebuffer_size_callback(GLFWwindow* window, int width, int height) override;
    virtual void mouse_move_callback(GLFWwindow* window, double xpos, double ypos) override;
    virtual void scroll_callback(GLFWwindow* window, double dx, double dy) override;
    virtual void
    mouse_button_callback(GLFWwindow* window, int button, int action, int mods) override;

    void render_objects(std::vector<common::renderable_node_t*> const& objects) const;

    void render_points();

    void transfer_vertices_to_gpu(
        unsigned int VBO,
        int position_attribute_location,
        int normal_attribute_location,
        int color_attribute_location,
        common::renderable_node_t const* object) const;

    void transfer_indices_to_gpu(unsigned int EBO, common::renderable_node_t const* object) const;

    void update_shader_view_projection_uniforms(shader_t const& shader) const;
    void update_shader_lighting_uniforms(shader_t const& shader) const;

    int get_position_attribute_location(shader_t const& shader) const;
    int get_normal_attribute_location(shader_t const& shader) const;
    int get_color_attribute_location(shader_t const& shader) const;

  private:
    void process_input(GLFWwindow* window, double dt);

  private:
    physics::xpbd::simulation_t& simulation_;
    point_light_t point_light_;
    directional_light_t directional_light_;

    camera_t camera_;

    GLFWwindow* window_;
    shader_t mesh_shader_;
    shader_t wireframe_shader_;
    shader_t point_shader_;

    std::vector<float> points_;
    bool should_render_points_;
    unsigned int point_vbo_;
    unsigned int point_vao_;

    std::vector<std::unique_ptr<common::renderable_node_t>> rendered_objects_;
};

} // namespace rendering
} // namespace sbs

/**
 * @brief Add helpful overloads for imgui inputs
 */
namespace ImGui {

template <typename Getter, typename Setter>
inline bool Checkbox(const char* label, Getter get, Setter set)
{
    bool value = get();
    bool ret   = ImGui::Checkbox(label, &value);
    set(value);
    return ret;
}

} // namespace ImGui

#endif // SBS_RENDERING_RENDERER_H