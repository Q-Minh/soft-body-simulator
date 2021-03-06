#include "sbs/rendering/renderer.h"

#include "sbs/rendering/pick.h"

#include <imgui/backends/imgui_impl_glfw.h>
#include <imgui/backends/imgui_impl_opengl3.h>
#include <imgui/imgui.h>
#include <sbs/physics/simulation.h>

namespace sbs {
namespace rendering {

renderer_base_t* renderer_base_t::active_renderer = nullptr;

renderer_t::renderer_t(
    physics::simulation_t& simulation,
    point_light_t const& point_light,
    directional_light_t const& directional_light)
    : simulation_(simulation),
      point_light_(point_light),
      directional_light_(directional_light),
      camera_(),
      window_(),
      mesh_shader_(),
      wireframe_shader_(),
      point_shader_(),
      points_(),
      should_render_points_(),
      point_vbo_(),
      point_vao_()
{
    initialize();
}

renderer_t::~renderer_t()
{
    close();
}

bool renderer_t::initialize()
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
    glfwSetMouseButtonCallback(window, renderer_base_t::mouse_button_callback_dispatcher);
    glfwSetCursorPosCallback(window, renderer_base_t::mouse_move_callback_dispatcher);
    glfwSetScrollCallback(window, renderer_base_t::scroll_callback_dispatcher);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        return false;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();

    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");

    window_ = window;

    glGenVertexArrays(1, &point_vao_);
    glGenBuffers(1, &point_vbo_);

    auto const create_objects_for_opengl =
        [](std::vector<std::unique_ptr<physics::body_t>> const& bodies) {
            for (std::unique_ptr<physics::body_t> const& body : bodies)
            {
                common::renderable_node_t& object = body->visual_model();
                unsigned int& VAO                 = object.VAO();
                unsigned int& VBO                 = object.VBO();
                unsigned int& EBO                 = object.EBO();

                glGenVertexArrays(1, &VAO);
                glGenBuffers(1, &VBO);
                glGenBuffers(1, &EBO);

                object.mark_vertices_dirty();
                object.mark_indices_dirty();
            }
        };

    create_objects_for_opengl(simulation_.bodies());

    return true;
}

bool renderer_t::use_mesh_shaders(
    std::filesystem::path const& vertex_shader_path,
    std::filesystem::path const& fragment_shader_path)
{
    mesh_shader_ = shader_t{vertex_shader_path, fragment_shader_path};
    return mesh_shader_.should_use();
}

bool renderer_t::use_wireframe_shaders(
    std::filesystem::path const& vertex_shader_path,
    std::filesystem::path const& fragment_shader_path)
{
    wireframe_shader_ = shader_t{vertex_shader_path, fragment_shader_path};
    return wireframe_shader_.should_use();
}

bool renderer_t::use_point_shaders(
    std::filesystem::path const& vertex_shader_path,
    std::filesystem::path const& fragment_shader_path)
{
    point_shader_ = shader_t(vertex_shader_path, fragment_shader_path);
    return point_shader_.should_use();
}

void renderer_t::launch()
{
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_POINT_SMOOTH);

    double last_frame_time = 0.f;
    while (!glfwWindowShouldClose(window_))
    {
        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        double const now = glfwGetTime();
        double const dt  = now - last_frame_time;
        last_frame_time  = now;

        process_input(window_, dt);

        if (on_new_imgui_frame)
        {
            on_new_imgui_frame(simulation_);
        }

        if (on_new_physics_timestep)
        {
            on_new_physics_timestep(dt, simulation_);
        }

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        if (on_pre_render)
        {
            on_pre_render(simulation_);
        }

        std::vector<common::renderable_node_t*> objects{};
        std::transform(
            simulation_.bodies().begin(),
            simulation_.bodies().end(),
            std::back_inserter(objects),
            [](std::unique_ptr<physics::body_t>& body) { return &body->visual_model(); });
        render_objects(objects);

        if (should_render_points_)
            render_points();

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window_);
    }
}

void renderer_t::close()
{
    auto const delete_objects_from_opengl =
        [](std::vector<std::unique_ptr<physics::body_t>> const& bodies) {
            for (std::unique_ptr<physics::body_t> const& body : bodies)
            {
                common::renderable_node_t const& object = body->visual_model();
                glDeleteVertexArrays(1, &(object.VAO()));
                glDeleteBuffers(1, &(object.VBO()));
                glDeleteBuffers(1, &(object.EBO()));
            }
        };

    delete_objects_from_opengl(simulation_.bodies());

    mesh_shader_.destroy();
    wireframe_shader_.destroy();
    point_shader_.destroy();
    glfwTerminate();
}

void renderer_t::add_point(std::array<float, 9u> const& xyz_nxnynz_rgb_point)
{
    std::copy(
        xyz_nxnynz_rgb_point.begin(),
        xyz_nxnynz_rgb_point.end(),
        std::back_inserter(points_));
    should_render_points_ = true;
}

void renderer_t::clear_points()
{
    points_.clear();
    should_render_points_ = true;
}

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

void renderer_t::mouse_move_callback(GLFWwindow* window, double xpos, double ypos)
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

    bool should_return{false};

    bool const is_picking = std::any_of(pickers.begin(), pickers.end(), [](picker_t const& picker) {
        return picker.is_picking();
    });
    should_return |= is_picking;

    if (is_picking && !ImGui::GetIO().WantCaptureMouse)
    {
        std::for_each(pickers.begin(), pickers.end(), [window, xpos, ypos](picker_t& picker) {
            if (picker.is_usable())
                picker.mouse_moved_event(window, xpos, ypos);
        });
    }
    if (on_mouse_moved && !ImGui::GetIO().WantCaptureMouse)
    {
        should_return |= on_mouse_moved(window, xpos, ypos);
    }
    if (should_return)
    {
        return;
    }

    int state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
    if (state == GLFW_PRESS && !ImGui::GetIO().WantCaptureMouse)
    {
        camera_.handle_mouse_movement(dx, dy);
    }
}

void renderer_t::scroll_callback(GLFWwindow* window, double dx, double dy)
{
    if (!ImGui::GetIO().WantCaptureMouse)
    {
        camera_.handle_mouse_scroll(static_cast<float>(dy));
    }
}

void renderer_t::mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (!ImGui::GetIO().WantCaptureMouse)
    {
        if (on_mouse_button_pressed)
        {
            on_mouse_button_pressed(window, button, action, mods);
        }
        std::for_each(
            pickers.begin(),
            pickers.end(),
            [window, button, action, mods](picker_t& picker) {
                if (picker.is_usable())
                    picker.mouse_button_pressed_event(window, button, action, mods);
            });
    }
}

void renderer_t::render_objects(std::vector<common::renderable_node_t*> const& objects) const
{
    /**
     * Render the scene. Transfers data to the GPU every frame, since
     * we are primarily working with soft bodies. This is of course
     * not optimal, but is great for prototyping.
     */
    for (common::renderable_node_t* object : objects)
    {
        bool const should_render_wireframe = object->should_render_wireframe();

        shader_t const& shader = should_render_wireframe ? wireframe_shader_ : mesh_shader_;
        shader.use();

        auto const position_attribute_location = get_position_attribute_location(shader);
        auto const normal_attribute_location   = get_normal_attribute_location(shader);
        auto const color_attribute_location    = get_color_attribute_location(shader);

        /**
         * Setup lights and view/projection
         */
        update_shader_view_projection_uniforms(shader);
        if (should_render_wireframe)
        {
            // wireframe shader does not need lighting uniforms as it performs flat shading
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        }
        else // should_render_triangles
        {
            update_shader_lighting_uniforms(shader);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        }

        unsigned int& VAO = object->VAO();
        unsigned int& VBO = object->VBO();
        unsigned int& EBO = object->EBO();

        glBindVertexArray(VAO);

        /**
         * Only send data to the GPU if geometry has changed
         */
        if (object->should_transfer_vertices())
        {
            object->prepare_vertices_for_rendering();
            transfer_vertices_to_gpu(
                VBO,
                position_attribute_location,
                normal_attribute_location,
                color_attribute_location,
                object);
            object->mark_vertices_clean();
        }

        if (object->should_transfer_indices())
        {
            object->prepare_indices_for_rendering();
            transfer_indices_to_gpu(EBO, object);
            object->mark_indices_clean();
        }

        auto const num_indices = object->get_cpu_index_buffer().size();

        /**
         * Draw the mesh
         */
        glDrawElements(GL_TRIANGLES, static_cast<int>(num_indices), GL_UNSIGNED_INT, 0);

        /**
         * Unbind buffers and vertex arrays
         */
        glBindBuffer(GL_ARRAY_BUFFER, 0u);
        glBindVertexArray(0u);
    }
}

void renderer_t::render_points()
{
    shader_t const& shader = point_shader_;
    shader.use();

    auto const position_attribute_location = get_position_attribute_location(shader);
    auto const normal_attribute_location   = get_normal_attribute_location(shader);
    auto const color_attribute_location    = get_color_attribute_location(shader);

    /**
     * Setup lights and view/projection
     */
    update_shader_view_projection_uniforms(shader);

    glBindVertexArray(point_vao_);
    glBindBuffer(GL_ARRAY_BUFFER, point_vbo_);

    auto constexpr num_bytes_per_float = sizeof(float);
    glBufferData(
        GL_ARRAY_BUFFER,
        points_.size() * num_bytes_per_float,
        points_.data(),
        GL_DYNAMIC_DRAW);

    auto const num_vertices           = points_.size() / 9u;
    auto constexpr size_of_one_vertex = 3u * num_bytes_per_float /* x,y,z coordinates */ +
                                        3u * num_bytes_per_float /* nx,ny,nz normal components */ +
                                        3u * num_bytes_per_float /* r,g,b colors */;
    auto constexpr stride_between_vertices = size_of_one_vertex;
    auto constexpr vertex_position_offset  = 0u;
    auto constexpr vertex_normal_offset    = vertex_position_offset + 3u * num_bytes_per_float;
    auto constexpr vertex_color_offset     = vertex_normal_offset + 3u * num_bytes_per_float;

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

    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(num_vertices));

    should_render_points_ = false;
}

void renderer_t::transfer_vertices_to_gpu(
    unsigned int VBO,
    int position_attribute_location,
    int normal_attribute_location,
    int color_attribute_location,
    common::renderable_node_t const* object) const
{
    auto constexpr num_bytes_per_float = sizeof(float);
    auto constexpr size_of_one_vertex  = 3u * num_bytes_per_float /* x,y,z coordinates */ +
                                        3u * num_bytes_per_float /* nx,ny,nz normal components */ +
                                        3u * num_bytes_per_float /* r,g,b colors */;

    std::vector<float> const& cpu_buffer = object->get_cpu_vertex_buffer();

    /**
     * Transfer vertex data to the GPU
     */
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(
        GL_ARRAY_BUFFER,
        cpu_buffer.size() * num_bytes_per_float,
        cpu_buffer.data(),
        GL_DYNAMIC_DRAW);

    /**
     * Specify vertex attributes' layout on the GPU
     */
    auto constexpr stride_between_vertices = size_of_one_vertex;
    auto constexpr vertex_position_offset  = 0u;
    auto constexpr vertex_normal_offset    = vertex_position_offset + 3u * num_bytes_per_float;
    auto constexpr vertex_color_offset     = vertex_normal_offset + 3u * num_bytes_per_float;

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
}

void renderer_t::transfer_indices_to_gpu(unsigned int EBO, common::renderable_node_t const* object)
    const
{
    std::vector<std::uint32_t> const& indices = object->get_cpu_index_buffer();
    auto const num_indices                    = indices.size();

    /**
     * Transfer triangle data to the GPU
     */
    auto constexpr num_bytes_per_index = sizeof(decltype(indices.front()));

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER,
        indices.size() * num_bytes_per_index,
        indices.data(),
        GL_DYNAMIC_DRAW);
}

void renderer_t::update_shader_view_projection_uniforms(shader_t const& shader) const
{
    int width  = 0;
    int height = 0;
    glfwGetWindowSize(window_, &width, &height);
    float const aspect_ratio   = static_cast<float>(width) / static_cast<float>(height);
    glm::mat4 const projection = camera_.projection_gl(aspect_ratio);
    glm::mat4 const view       = camera_.view_gl();

    shader.set_mat4_uniform("projection", projection);
    shader.set_mat4_uniform("view", view);
}

void renderer_t::update_shader_lighting_uniforms(shader_t const& shader) const
{
    auto const& directional_light = directional_light_;
    auto const& point_light       = point_light_;

    shader.set_vec3_uniform("ViewPosition", camera_.position());

    shader.set_vec3_uniform(
        "DirectionalLight.direction",
        glm::vec3{directional_light.dx, directional_light.dy, directional_light.dz});
    shader.set_vec3_uniform(
        "DirectionalLight.ambient",
        glm::vec3{
            directional_light.ambient.r,
            directional_light.ambient.g,
            directional_light.ambient.b});
    shader.set_vec3_uniform(
        "DirectionalLight.diffuse",
        glm::vec3{
            directional_light.diffuse.r,
            directional_light.diffuse.g,
            directional_light.diffuse.b});
    shader.set_vec3_uniform(
        "DirectionalLight.specular",
        glm::vec3{
            directional_light.specular.r,
            directional_light.specular.g,
            directional_light.specular.b});
    shader.set_float_uniform("DirectionalLight.exponent", directional_light.specular.exp);

    shader.set_vec3_uniform(
        "PointLight.position",
        glm::vec3{point_light.x, point_light.y, point_light.z});
    shader.set_vec3_uniform(
        "PointLight.ambient",
        glm::vec3{point_light.ambient.r, point_light.ambient.g, point_light.ambient.b});
    shader.set_vec3_uniform(
        "PointLight.diffuse",
        glm::vec3{point_light.diffuse.r, point_light.diffuse.g, point_light.diffuse.b});
    shader.set_vec3_uniform(
        "PointLight.specular",
        glm::vec3{point_light.specular.r, point_light.specular.g, point_light.specular.b});
    shader.set_float_uniform("PointLight.exponent", point_light.specular.exp);
    shader.set_float_uniform("PointLight.constant", point_light.attenuation.constant);
    shader.set_float_uniform("PointLight.linear", point_light.attenuation.linear);
    shader.set_float_uniform("PointLight.quadratic", point_light.attenuation.quadratic);
}

int renderer_t::get_position_attribute_location(shader_t const& shader) const
{
    return glGetAttribLocation(shader.id(), shader_t::vertex_shader_position_attribute_name);
}

int renderer_t::get_normal_attribute_location(shader_t const& shader) const
{
    return glGetAttribLocation(shader.id(), shader_t::vertex_shader_normal_attribute_name);
}

int renderer_t::get_color_attribute_location(shader_t const& shader) const
{
    return glGetAttribLocation(shader.id(), shader_t::vertex_shader_color_attribute_name);
}

} // namespace rendering
} // namespace sbs
