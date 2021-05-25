#ifndef SBS_RENDERING_RENDERER_H
#define SBS_RENDERING_RENDERER_H

// clang-format off
#include <glad/glad.h>
#include <GLFW/glfw3.h>
// clang-format on

#include "io/load_scene.h"
#include "camera.h"
#include "shader.h"

#include <filesystem>
#include <functional>
#include <iostream>
#include <iterator>

namespace sbs {
namespace rendering {

class renderer_base_t
{
  public:
    virtual void framebuffer_size_callback(GLFWwindow* window, int width, int height) = 0;
    virtual void mouse_callback(GLFWwindow* window, double xpos, double ypos)         = 0;
    virtual void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)  = 0;

    static renderer_base_t* active_renderer;

    virtual void set_as_active_renderer() { active_renderer = this; }

    static void framebuffer_size_callback_dispatcher(GLFWwindow* window, int width, int height)
    {
        if (active_renderer != nullptr)
            active_renderer->framebuffer_size_callback(window, width, height);
    }
    static void mouse_callback_dispatcher(GLFWwindow* window, double xpos, double ypos)
    {
        if (active_renderer != nullptr)
            active_renderer->mouse_callback(window, xpos, ypos);
    }
    static void scroll_callback_dispatcher(GLFWwindow* window, double xoffset, double yoffset)
    {
        if (active_renderer != nullptr)
            active_renderer->scroll_callback(window, xoffset, yoffset);
    }
};

class renderer_t : public renderer_base_t
{
  public:
    bool initialize(std::filesystem::path const& scene_path);

    bool use_shaders(
        std::filesystem::path const& vertex_shader_path,
        std::filesystem::path const& fragment_shader_path);

    void launch();

    void close();

    std::vector<std::string> get_error_messages() const { return shader_.error_messages(); }

    virtual void framebuffer_size_callback(GLFWwindow* window, int width, int height) override;
    virtual void mouse_callback(GLFWwindow* window, double xpos, double ypos) override;
    virtual void scroll_callback(GLFWwindow* window, double dx, double dy) override;

    std::uint32_t constexpr get_initial_window_width() const { return 800u; }
    std::uint32_t constexpr get_initial_window_height() const { return 600u; }

    std::function<void(common::scene_t&)> on_scene_loaded;
    std::function<void(double /*render_frame_dt*/, common::scene_t& /*scene_*/)> on_new_frame;

  private:
    void process_input(GLFWwindow* window, double dt);

  private:
    common::scene_t scene_;
    camera_t camera_;

    GLFWwindow* window_;
    shader_t shader_;
};

} // namespace rendering
} // namespace sbs

#endif // SBS_RENDERING_RENDERER_H