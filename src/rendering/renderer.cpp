#include "rendering/renderer.h"

#include <iostream>

namespace sbs {
namespace rendering {

renderer_base_t* renderer_base_t::active_renderer = nullptr;

void renderer_t::process_input(GLFWwindow* window, double dt)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, true);
    }

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
    {
        if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
        {
            camera_.handle_keyboard(camera_t::movement_t::up, static_cast<float>(dt));
        }
        else
        {
            camera_.handle_keyboard(camera_t::movement_t::forward, static_cast<float>(dt));
        }
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
    {
        camera_.handle_keyboard(camera_t::movement_t::left, static_cast<float>(dt));
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
    {
        if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS)
        {
            camera_.handle_keyboard(camera_t::movement_t::down, static_cast<float>(dt));
        }
        else
        {
            camera_.handle_keyboard(camera_t::movement_t::backward, static_cast<float>(dt));
        }
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
    {
        camera_.handle_keyboard(camera_t::movement_t::right, static_cast<float>(dt));
    }
}

void renderer_t::framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void renderer_t::mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    static bool first_mouse_movement = true;
    static double last_x_pos         = get_initial_window_width() / 2.f;
    static double last_y_pos         = get_initial_window_height() / 2.f;

    if (first_mouse_movement)
    {
        last_x_pos           = xpos;
        last_y_pos           = ypos;
        first_mouse_movement = false;
    }

    float const dx = static_cast<float>(xpos - last_x_pos);
    float const dy = static_cast<float>(last_y_pos - ypos);

    last_x_pos = xpos;
    last_y_pos = ypos;

    int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);

    if (state == GLFW_PRESS)
    {
        camera_.handle_mouse_movement(dx, dy);
    }
}

void renderer_t::scroll_callback(GLFWwindow* window, double dx, double dy)
{
    camera_.handle_mouse_scroll(static_cast<float>(dy));
}

bool renderer_t::initialize(std::filesystem::path const& scene_path)
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    int const width    = static_cast<int>(get_initial_window_width());
    int const height   = static_cast<int>(get_initial_window_height());
    GLFWwindow* window = glfwCreateWindow(width, height, "Soft Body Simulator", NULL, NULL);
    if (window == NULL)
    {
        return false;
    }
    glfwMakeContextCurrent(window);

    this->set_as_active_renderer();
    glfwSetFramebufferSizeCallback(window, renderer_base_t::framebuffer_size_callback_dispatcher);
    // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetCursorPosCallback(window, renderer_base_t::mouse_callback_dispatcher);
    glfwSetScrollCallback(window, renderer_base_t::scroll_callback_dispatcher);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        return false;
    }

    window_ = window;

    scene_ = io::load_scene(scene_path);
    for (auto& object : scene_.objects)
    {
        unsigned int& VAO = object->VAO;
        unsigned int& VBO = object->VBO;
        unsigned int& EBO = object->EBO;

        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glGenBuffers(1, &EBO);
    }

    if (on_scene_loaded)
    {
        on_scene_loaded(scene_);
    }

    return true;
}

bool renderer_t::use_shaders(
    std::filesystem::path const& vertex_shader_path,
    std::filesystem::path const& fragment_shader_path)
{
    shader_ = shader_t{vertex_shader_path, fragment_shader_path};
    return shader_.should_use();
}

void renderer_t::launch()
{
    glEnable(GL_DEPTH_TEST);

    auto const position_attribute_location =
        glGetAttribLocation(shader_.id(), shader_t::vertex_shader_position_attribute_name);
    auto const normal_attribute_location =
        glGetAttribLocation(shader_.id(), shader_t::vertex_shader_normal_attribute_name);
    auto const color_attribute_location =
        glGetAttribLocation(shader_.id(), shader_t::vertex_shader_color_attribute_name);

    double last_frame_time = 0.f;
    while (!glfwWindowShouldClose(window_))
    {
        double const now = glfwGetTime();
        double const dt  = now - last_frame_time;
        last_frame_time  = now;

        process_input(window_, dt);

        if (on_new_imgui_frame)
        {
            on_new_imgui_frame(scene_);
        }

        if (on_new_physics_timestep)
        {
            on_new_physics_timestep(dt, scene_);
        }

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        shader_.use();

        int width  = 0;
        int height = 0;
        glfwGetWindowSize(window_, &width, &height);
        float const aspect_ratio   = static_cast<float>(width) / static_cast<float>(height);
        glm::mat4 const projection = camera_.projection_matrix(aspect_ratio);
        glm::mat4 const view       = camera_.view_matrix();

        shader_.set_mat4_uniform("projection", projection);
        shader_.set_mat4_uniform("view", view);

        auto const& directional_light = scene_.directional_light;
        auto const& point_light       = scene_.point_light;

        shader_.set_vec3_uniform("ViewPosition", camera_.position());

        shader_.set_vec3_uniform(
            "DirectionalLight.direction",
            glm::vec3{directional_light.dx, directional_light.dy, directional_light.dz});
        shader_.set_vec3_uniform(
            "DirectionalLight.ambient",
            glm::vec3{
                directional_light.ambient.r,
                directional_light.ambient.g,
                directional_light.ambient.b});
        shader_.set_vec3_uniform(
            "DirectionalLight.diffuse",
            glm::vec3{
                directional_light.diffuse.r,
                directional_light.diffuse.g,
                directional_light.diffuse.b});
        shader_.set_vec3_uniform(
            "DirectionalLight.specular",
            glm::vec3{
                directional_light.specular.r,
                directional_light.specular.g,
                directional_light.specular.b});
        shader_.set_float_uniform("DirectionalLight.exponent", directional_light.specular.exp);

        shader_.set_vec3_uniform(
            "PointLight.position",
            glm::vec3{point_light.x, point_light.y, point_light.z});
        shader_.set_vec3_uniform(
            "PointLight.ambient",
            glm::vec3{point_light.ambient.r, point_light.ambient.g, point_light.ambient.b});
        shader_.set_vec3_uniform(
            "PointLight.diffuse",
            glm::vec3{point_light.diffuse.r, point_light.diffuse.g, point_light.diffuse.b});
        shader_.set_vec3_uniform(
            "PointLight.specular",
            glm::vec3{point_light.specular.r, point_light.specular.g, point_light.specular.b});
        shader_.set_float_uniform("PointLight.exponent", point_light.specular.exp);
        shader_.set_float_uniform("PointLight.constant", point_light.attenuation.constant);
        shader_.set_float_uniform("PointLight.linear", point_light.attenuation.linear);
        shader_.set_float_uniform("PointLight.quadratic", point_light.attenuation.quadratic);

        /**
         * Render the scene. Transfers data to the GPU every frame, since
         * we are primarily working with soft bodies. This is of course
         * not optimal, but is great for prototyping.
         */
        for (auto& object : scene_.objects)
        {
            unsigned int& VAO = object->VAO;
            unsigned int& VBO = object->VBO;
            unsigned int& EBO = object->EBO;

            object->mesh.extract_boundary_surface_mesh();

            glBindVertexArray(VAO);

            std::vector<float> cpu_buffer{};
            auto constexpr num_bytes_per_float = sizeof(float);
            auto constexpr size_of_one_vertex =
                3u * num_bytes_per_float /* x,y,z coordinates */ +
                3u * num_bytes_per_float /* nx,ny,nz normal components */ +
                3u * num_bytes_per_float /* r,g,b colors */;

            auto const number_of_vertices =
                static_cast<std::size_t>(object->mesh.vertices().cols());
            cpu_buffer.reserve(number_of_vertices * size_of_one_vertex);
            for (std::size_t i = 0u; i < number_of_vertices; ++i)
            {
                float const x = static_cast<float>(object->mesh.vertices()(0u, i));
                float const y = static_cast<float>(object->mesh.vertices()(1u, i));
                float const z = static_cast<float>(object->mesh.vertices()(2u, i));
                cpu_buffer.push_back(x);
                cpu_buffer.push_back(y);
                cpu_buffer.push_back(z);

                float const nx = static_cast<float>(object->mesh.normals()(0u, i));
                float const ny = static_cast<float>(object->mesh.normals()(1u, i));
                float const nz = static_cast<float>(object->mesh.normals()(2u, i));
                cpu_buffer.push_back(nx);
                cpu_buffer.push_back(ny);
                cpu_buffer.push_back(nz);

                float const r = static_cast<float>(object->mesh.colors()(0u, i));
                float const g = static_cast<float>(object->mesh.colors()(1u, i));
                float const b = static_cast<float>(object->mesh.colors()(2u, i));
                cpu_buffer.push_back(r);
                cpu_buffer.push_back(g);
                cpu_buffer.push_back(b);
            }

            auto const number_of_faces = static_cast<std::size_t>(object->mesh.faces().cols());
            std::vector<std::uint32_t> indices{};
            indices.reserve(number_of_faces);
            for (std::size_t f = 0u; f < number_of_faces; ++f)
            {
                indices.push_back(object->mesh.faces()(0u, f));
                indices.push_back(object->mesh.faces()(1u, f));
                indices.push_back(object->mesh.faces()(2u, f));
            }

            /**
             * Transfer data to the GPU
             */
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(
                GL_ARRAY_BUFFER,
                cpu_buffer.size() * num_bytes_per_float,
                cpu_buffer.data(),
                GL_DYNAMIC_DRAW);

            auto constexpr num_bytes_per_index = sizeof(decltype(indices.front()));

            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
            glBufferData(
                GL_ELEMENT_ARRAY_BUFFER,
                indices.size() * num_bytes_per_index,
                indices.data(),
                GL_DYNAMIC_DRAW);

            /**
             * Specify vertex attributes' layout on the GPU
             */

            auto constexpr stride_between_vertices = size_of_one_vertex;
            auto constexpr vertex_position_offset  = 0u;
            auto constexpr vertex_normal_offset    = 3u * num_bytes_per_float;
            auto constexpr vertex_color_offset = vertex_normal_offset + 3u * num_bytes_per_float;

            glVertexAttribPointer(
                position_attribute_location,
                3u,
                GL_FLOAT,
                GL_FALSE,
                stride_between_vertices,
                reinterpret_cast<void*>(vertex_position_offset));
            glEnableVertexAttribArray(position_attribute_location);

            glVertexAttribPointer(
                normal_attribute_location,
                3u,
                GL_FLOAT,
                GL_FALSE,
                stride_between_vertices,
                reinterpret_cast<void*>(vertex_normal_offset));
            glEnableVertexAttribArray(normal_attribute_location);

            glVertexAttribPointer(
                color_attribute_location,
                3u,
                GL_FLOAT,
                GL_FALSE,
                stride_between_vertices,
                reinterpret_cast<void*>(vertex_color_offset));
            glEnableVertexAttribArray(color_attribute_location);

            /**
             * Draw the mesh
             */
            glDrawElements(GL_TRIANGLES, static_cast<int>(indices.size()), GL_UNSIGNED_INT, 0);

            /**
             * Unbind buffers and vertex arrays
             */
            glBindBuffer(GL_ARRAY_BUFFER, 0u);
            glBindVertexArray(0u);
        }

        glfwSwapBuffers(window_);
        glfwPollEvents();
    }
}

void renderer_t::close()
{
    for (auto const& object : scene_.objects)
    {
        glDeleteVertexArrays(1, &(object->VAO));
        glDeleteBuffers(1, &(object->VBO));
        glDeleteBuffers(1, &(object->EBO));
    }
    shader_.destroy();
    glfwTerminate();
}

} // namespace rendering
} // namespace sbs
