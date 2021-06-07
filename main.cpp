#include "io/load_scene.h"
#include "physics/cutting/cut_tetrahedron.h"
#include "physics/xpbd/solver.h"
#include "rendering/pick.h"
#include "rendering/renderer.h"
#include "rendering/trackball_rotation_adapter.h"

#include <chrono>
#include <imgui/imgui.h>
#include <iostream>

int main(int argc, char** argv)
{
    std::filesystem::path const cwd = std::filesystem::current_path();

    if (argc != 4)
    {
        std::cerr << "Usage: sbs-viewer.exe <scene specification json file> "
                     "<path/to/vertex_shader.vs> <path/to/fragment_shader.fs>\n";
        return 1;
    }

    std::filesystem::path const scene_specification_path{argv[1]};
    std::filesystem::path const vertex_shader_path{argv[2]};
    std::filesystem::path const fragment_shader_path{argv[3]};

    sbs::rendering::renderer_t renderer{};
    sbs::physics::xpbd::solver_t solver{};

    auto const initial_scene = sbs::io::load_scene(scene_specification_path);

    auto const cutting_surface_node_it = std::find_if(
        initial_scene.environment_objects.begin(),
        initial_scene.environment_objects.end(),
        [](std::shared_ptr<sbs::common::node_t> const object) {
            return object->id == "cutting surface";
        });

    sbs::rendering::trackball_rotation_adapter_t cutting_surface_trackball_adapter{};
    if (cutting_surface_node_it != initial_scene.environment_objects.end())
    {
        std::shared_ptr<sbs::common::node_t> const cutting_surface_node = *cutting_surface_node_it;
        sbs::common::shared_vertex_surface_mesh_t model_space_cutting_surface =
            cutting_surface_node->render_model;

        auto const edges = model_space_cutting_surface.boundary_edges();

        cutting_surface_trackball_adapter =
            sbs::rendering::trackball_rotation_adapter_t{model_space_cutting_surface};
    }

    /**
     * physics update goes here
     */
    std::vector<sbs::physics::xpbd::simulation_parameters_t> per_body_simulation_parameters{};
    renderer.on_scene_loaded = [&](sbs::common::scene_t& scene) {
        per_body_simulation_parameters.clear();
        for (auto const& body : scene.physics_objects)
        {
            sbs::physics::xpbd::simulation_parameters_t params{};
            per_body_simulation_parameters.push_back(params);

            body->physical_model.forces().setZero();
        }

        solver.setup(&scene.physics_objects, per_body_simulation_parameters);
    };

    bool are_physics_active          = false;
    double timestep                  = 1. / 60.;
    std::uint32_t iterations         = 60u;
    std::uint32_t substeps           = 60u;
    double fps                       = 0u;
    renderer.on_new_physics_timestep = [&](double render_frame_dt, sbs::common::scene_t& scene) {
        static double tb = 0.;

        tb += render_frame_dt;
        auto const time_between_frames = tb;

        /**
         * If the elapsed time between the last physics update was less than
         * the physics timestep, we don't update physics.
         */
        if (tb < timestep)
            return;

        double const num_timesteps_elapsed = std::floor(tb / timestep);
        tb -= num_timesteps_elapsed * timestep;

        auto const begin = std::chrono::steady_clock::now();

        /**
         * Compute external forces
         */
        if (are_physics_active)
        {
            for (auto const& body : scene.physics_objects)
            {
                /**
                 * Reset external forces
                 */
                body->physical_model.forces().colwise() += Eigen::Vector3d{0., -9.81, 0.};
            }

            solver.step(timestep, iterations, substeps);

            for (auto const& body : scene.physics_objects)
            {
                body->render_model = body->physical_model.boundary_surface_mesh().to_face_based();
                body->render_state.should_transfer_vertices = true;
                body->render_state.should_transfer_indices  = true;
                body->physical_model.forces().setZero();
            }
        }

        auto const end = std::chrono::steady_clock::now();
        auto const duration =
            std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();

        fps = 1. / static_cast<double>(duration) * 1e9;
    };

    int active_environment_body_idx = 0;
    int active_physics_body_idx     = 0;

    renderer.on_mouse_button_pressed = [&](GLFWwindow* window, int button, int action, int mods) {
        bool const is_picking =
            (button == GLFW_MOUSE_BUTTON_LEFT && mods == GLFW_MOD_SHIFT && action == GLFW_PRESS);
        if (!is_picking)
            return;

        sbs::common::scene_t const& scene = renderer.scene();
        sbs::common::node_t& active_node  = *scene.physics_objects[active_physics_body_idx];
        sbs::common::shared_vertex_surface_mesh_t const& surface =
            active_node.physical_model.boundary_surface_mesh();

        int viewport_gl[4];
        glGetIntegerv(GL_VIEWPORT, viewport_gl);
        Eigen::Vector4d const viewport{
            static_cast<double>(viewport_gl[0]),
            static_cast<double>(viewport_gl[1]),
            static_cast<double>(viewport_gl[2]),
            static_cast<double>(viewport_gl[3])};

        float const aspect_ratio =
            static_cast<float>(viewport(2)) / static_cast<float>(viewport(3));
        Eigen::Matrix4d const projection = renderer.camera().projection_matrix(aspect_ratio);
        Eigen::Matrix4d const view       = renderer.camera().view_matrix();

        double x, y;
        glfwGetCursorPos(window, &x, &y);

        auto const ray = sbs::rendering::unproject_ray({x, y}, viewport, projection, view);
        auto const vi  = sbs::rendering::pick_vertex(ray, surface);
        if (vi.has_value())
        {
            active_node.physical_model.masses()(*vi) = 1e15;
        }
    };

    double xprev = 0., yprev = 0., dx = 0., dy = 0.;
    bool is_first_mouse_move_callback           = true;
    bool should_handle_cutting_surface_rotation = false;
    renderer.on_mouse_moved = [&](GLFWwindow* window, double x, double y) -> bool {
        bool result = false;

        int left_ctrl_key_state     = glfwGetKey(window, GLFW_KEY_LEFT_CONTROL);
        int left_mouse_button_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);

        if (left_ctrl_key_state == GLFW_PRESS && left_mouse_button_state == GLFW_PRESS)
        {
            dx                                     = x - xprev;
            dy                                     = yprev - y;
            should_handle_cutting_surface_rotation = true;
            result                                 = true;
        }

        xprev = x;
        yprev = y;

        return result;
    };

    renderer.on_new_imgui_frame = [&](sbs::common::scene_t& scene) {
        ImGui::Begin("Soft Body Simulator");

        std::shared_ptr<sbs::common::node_t> active_environment_node =
            scene.environment_objects[active_environment_body_idx];

        std::shared_ptr<sbs::common::node_t> active_physics_node =
            scene.physics_objects[active_physics_body_idx];

        if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_None))
        {
            ImGui::TreePush();
            float const scene_panel_width = ImGui::GetColumnWidth();
            if (ImGui::CollapsingHeader("Environment Objects##Scene", ImGuiTreeNodeFlags_None))
            {
                ImGui::BulletText("Select active object");
                for (std::size_t b = 0u; b < scene.environment_objects.size(); ++b)
                {
                    auto const& body = scene.environment_objects[b];
                    ImGui::RadioButton(
                        body->id.c_str(),
                        &active_environment_body_idx,
                        static_cast<int>(b));
                }

                active_environment_node = scene.environment_objects[active_environment_body_idx];
            }
            if (ImGui::CollapsingHeader("Physics Objects##Scene", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ImGui::BulletText("Select active object");
                for (std::size_t b = 0u; b < scene.physics_objects.size(); ++b)
                {
                    auto const& body = scene.physics_objects[b];
                    ImGui::RadioButton(
                        body->id.c_str(),
                        &active_physics_body_idx,
                        static_cast<int>(b));
                }

                active_physics_node = scene.physics_objects[active_physics_body_idx];
            }
            ImGui::TreePop();

            static bool is_wireframe = false;
            ImGui::Checkbox(
                "Wireframe",
                [&]() { return is_wireframe; },
                [&](bool want_render_wireframe) {
                    bool dirty = false;

                    if (!is_wireframe && want_render_wireframe)
                    {
                        for (auto const& body : scene.physics_objects)
                        {
                            body->render_model = body->physical_model.facets().to_face_based();
                            body->render_state.should_render_wireframe = true;
                        }
                        dirty = true;
                    }
                    if (is_wireframe && !want_render_wireframe)
                    {
                        for (auto const& body : scene.physics_objects)
                        {
                            body->render_model =
                                body->physical_model.boundary_surface_mesh().to_face_based();
                            body->render_state.should_render_wireframe = false;
                        }
                        dirty = true;
                    }

                    if (dirty)
                    {
                        for (auto const& body : scene.physics_objects)
                        {
                            body->render_state.should_transfer_vertices = true;
                            body->render_state.should_transfer_indices  = true;
                        }
                    }

                    is_wireframe = want_render_wireframe;
                });

            if (ImGui::Button("Reload", ImVec2(scene_panel_width / 2.f, 0.f)))
            {
                renderer.unload_current_scene();
                scene = sbs::io::load_scene(scene_specification_path);
                renderer.load_scene(scene);
                active_physics_body_idx     = 0;
                active_physics_node         = scene.physics_objects[active_physics_body_idx];
                active_environment_body_idx = 0;
                active_environment_node = scene.environment_objects[active_environment_body_idx];
            }
        }

        if (ImGui::CollapsingHeader("Physics", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::TreePush();
            if (ImGui::CollapsingHeader("XPBD"))
            {
                ImGui::Text("Constraint type");
                ImGui::TreePush();
                static int constraint_type_choice = 0;
                ImGui::RadioButton("Green##XPBD", &constraint_type_choice, 0);
                ImGui::RadioButton("Distance##XPBD", &constraint_type_choice, 1);
                ImGui::TreePop();

                static float alpha = 1e-8f;
                ImGui::InputFloat("Compliance##XPBD", &alpha, 0.0000001f, 0.1f, "%.8f");
                static float mass_per_vertex = 1.f;
                ImGui::InputFloat("Vertex mass##XPBD", &mass_per_vertex, 1.f, 10.f, "%.1f");

                ImGui::Text("FEM");
                static float young_modulus = 1e12f;
                static float poisson_ratio = 0.45f;
                ImGui::TreePush();
                ImGui::InputFloat("Young modulus##XPBD", &young_modulus, 1000.f, 10'000.f, "%.1f");
                ImGui::InputFloat("Poisson ratio##XPBD", &poisson_ratio, 0.01f, 0.1f, "%.2f");
                ImGui::TreePop();

                ImGui::Text("Solver##XPBD");
                ImGui::TreePush();
                static float _timestep = timestep;
                static int _iterations = iterations;
                static int _substeps   = substeps;
                ImGui::InputFloat("Timestep##XPBD", &_timestep, 0.001f, 0.01f);
                ImGui::InputInt("Iterations##XPBD", &_iterations);
                ImGui::InputInt("Substeps##XPBD", &_substeps);
                ImGui::TreePop();

                float const w = ImGui::GetColumnWidth();
                if (ImGui::Button("Apply##XPBD", ImVec2(w / 2.f, 0.f)))
                {
                    per_body_simulation_parameters[active_physics_body_idx].alpha =
                        static_cast<double>(alpha);
                    per_body_simulation_parameters[active_physics_body_idx].constraint_type =
                        constraint_type_choice == 0 ?
                            sbs::physics::xpbd::constraint_type_t::green :
                            sbs::physics::xpbd::constraint_type_t::distance;
                    per_body_simulation_parameters[active_physics_body_idx].young_modulus =
                        static_cast<double>(young_modulus);
                    per_body_simulation_parameters[active_physics_body_idx].poisson_ratio =
                        static_cast<double>(poisson_ratio);

                    timestep   = static_cast<double>(_timestep);
                    iterations = static_cast<std::uint32_t>(_iterations);
                    substeps   = static_cast<std::uint32_t>(_substeps);

                    solver.setup(&scene.physics_objects, per_body_simulation_parameters);
                }

                ImGui::Checkbox("Activate physics", &are_physics_active);
            }
            if (ImGui::CollapsingHeader("Mesh"))
            {
                std::size_t const num_vertices = static_cast<std::size_t>(
                    active_physics_node->physical_model.positions().cols());
                ImGui::Text("Vertices: %d", num_vertices);
                std::size_t const num_elements =
                    static_cast<std::size_t>(active_physics_node->physical_model.elements().cols());
                ImGui::Text("Tetrahedra: %d", num_elements);
            }
            ImGui::TreePop();
            std::string const fps_str = "FPS: " + std::to_string(fps);
            ImGui::Text(fps_str.c_str());
        }

        if (ImGui::CollapsingHeader("Cutting", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (active_environment_node->id == "cutting surface")
            {
                ImGui::TreePush();

                static float tx = 0.f, ty = 0.f, tz = 0.f;
                float const w = ImGui::GetColumnWidth();

                ImGui::Text("Cutting surface translation");
                ImGui::SetNextItemWidth(0.95f * w);
                ImGui::SliderFloat("x##Cutting", &tx, -20.f, 20.f, "%.2f");
                ImGui::SetNextItemWidth(0.95f * w);
                ImGui::SliderFloat("y##Cutting", &ty, -20.f, 20.f, "%.2f");
                ImGui::SetNextItemWidth(0.95f * w);
                ImGui::SliderFloat("z##Cutting", &tz, -20.f, 20.f, "%.2f");

                static float sensitivity  = 0.005f;
                static double yaw_angle   = 0.;
                static double pitch_angle = 0.;

                ImGui::Text("Cutting surface rotation");
                ImGui::SetNextItemWidth(0.6f * w);
                ImGui::SliderFloat("sensitivity##Cutting", &sensitivity, 0.000001f, 0.01f, "%.6f");

                Eigen::Vector3d const translation{
                    static_cast<double>(tx),
                    static_cast<double>(ty),
                    static_cast<double>(tz)};

                if (should_handle_cutting_surface_rotation)
                {
                    cutting_surface_trackball_adapter.set_rotation_speed(
                        static_cast<double>(sensitivity));

                    cutting_surface_trackball_adapter.set_yaw_axis(Eigen::Vector3d{
                        renderer.camera().front().x,
                        renderer.camera().front().y,
                        renderer.camera().front().z});

                    cutting_surface_trackball_adapter.set_pitch_axis(Eigen::Vector3d{
                        -renderer.camera().right().x,
                        -renderer.camera().right().y,
                        -renderer.camera().right().z});

                    cutting_surface_trackball_adapter.rotate(dx, dy);

                    should_handle_cutting_surface_rotation = false;
                }

                active_environment_node->render_model.vertices() =
                    cutting_surface_trackball_adapter.mesh().vertices();
                active_environment_node->render_model.vertices().colwise() += translation;
                active_environment_node->render_state.should_transfer_vertices = true;

                if (ImGui::Button("Cut##Cutting", ImVec2(w / 2.f, 0.f)))
                {
                    bool const is_tetrahedral_mesh =
                        (active_physics_node->physical_model.elements().rows() == 4);

                    if (is_tetrahedral_mesh)
                    {
                        bool const has_mesh_been_cut = sbs::physics::cutting::cut_tetrahedral_mesh(
                            active_physics_node->physical_model,
                            active_environment_node->render_model);

                        if (has_mesh_been_cut)
                        {
                            solver.notify_topology_changed();
                            active_physics_node->render_model =
                                active_physics_node->render_state.should_render_wireframe ?
                                    active_physics_node->physical_model.facets().to_face_based() :
                                    active_physics_node->physical_model.boundary_surface_mesh()
                                        .to_face_based();
                            active_physics_node->render_state.should_transfer_vertices = true;
                            active_physics_node->render_state.should_transfer_indices  = true;
                        }
                    }
                }

                ImGui::TreePop();
            }
        }

        ImGui::End();
    };

    bool const initialization_success = renderer.initialize();
    bool const shader_loading_success =
        renderer.use_shaders(vertex_shader_path, fragment_shader_path);

    if (initialization_success && shader_loading_success)
    {
        renderer.load_scene(initial_scene);
        renderer.launch();
    }
    else
    {
        for (auto const& error_message : renderer.get_error_messages())
        {
            std::cout << error_message << "\n";
        }
    }

    return 0;
}