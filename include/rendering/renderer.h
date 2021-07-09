#ifndef SBS_RENDERING_RENDERER_H
#define SBS_RENDERING_RENDERER_H

// clang-format off
#include <glad/glad.h>
#include <GLFW/glfw3.h>
// clang-format on

#include "common/scene.h"
#include "camera.h"
#include "shader.h"

#include <filesystem>
#include <functional>

namespace sbs {
namespace rendering {

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
    bool initialize();
    void unload_current_scene();
    void load_scene(common::scene_t const& scene);

    void remove_object_from_scene(std::uint32_t object_idx);
    std::uint32_t add_object_to_scene(std::shared_ptr<common::renderable_node_t> const& node);

    bool use_shaders(
        std::filesystem::path const& vertex_shader_path,
        std::filesystem::path const& fragment_shader_path);

    bool use_wireframe_shaders(
        std::filesystem::path const& vertex_shader_path,
        std::filesystem::path const& fragment_shader_path);

    void launch();
    void close();

    std::vector<std::string> get_error_messages() const { return shader_.error_messages(); }

    std::uint32_t constexpr get_initial_window_width() const { return 800u; }
    std::uint32_t constexpr get_initial_window_height() const { return 600u; }

    camera_t const& camera() const { return camera_; }
    common::scene_t const& scene() const { return scene_; }

    /**
     * Render loop hooks
     */
    std::function<void(common::scene_t&)> on_scene_loaded;
    std::function<void(common::scene_t&)> on_new_imgui_frame;
    std::function<void(double /*render_frame_dt*/, common::scene_t& /*scene*/)>
        on_new_physics_timestep;

    /**
     * User input hooks
     */
    std::function<bool(GLFWwindow* /*window*/, double /*mouse_x_pos*/, double /*mouse_y_pos*/)>
        on_mouse_moved;
    std::function<void(GLFWwindow*, int, int, int)> on_mouse_button_pressed;

  protected:
    virtual void framebuffer_size_callback(GLFWwindow* window, int width, int height) override;
    virtual void mouse_move_callback(GLFWwindow* window, double xpos, double ypos) override;
    virtual void scroll_callback(GLFWwindow* window, double dx, double dy) override;
    virtual void
    mouse_button_callback(GLFWwindow* window, int button, int action, int mods) override;

    void
    render_objects(std::vector<std::shared_ptr<common::renderable_node_t>> const& objects) const;

    void transfer_vertices_to_gpu(
        unsigned int VBO,
        int position_attribute_location,
        int normal_attribute_location,
        int color_attribute_location,
        std::shared_ptr<common::renderable_node_t> const& object) const;

    void transfer_indices_to_gpu(
        unsigned int EBO,
        std::shared_ptr<common::renderable_node_t> const& object) const;

    void update_shader_view_projection_uniforms(shader_t const& shader) const;
    void update_shader_lighting_uniforms(shader_t const& shader) const;

  private:
    void process_input(GLFWwindow* window, double dt);

  private:
    common::scene_t scene_;
    camera_t camera_;

    GLFWwindow* window_;
    shader_t shader_;
    shader_t wireframe_shader_;

    bool should_render_wireframe_ = false;
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