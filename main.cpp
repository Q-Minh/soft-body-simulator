#include "sbs/common/mesh.h"
#include "sbs/io/load_scene.h"
#include "sbs/physics/xpbd/mesh.h"
#include "sbs/physics/xpbd/solver.h"
#include "sbs/rendering/physics_timestep_throttler.h"
#include "sbs/rendering/pick.h"
#include "sbs/rendering/renderer.h"

#include <imgui/imgui.h>
#include <iostream>

int main(int argc, char** argv)
{
    std::filesystem::path const cwd = std::filesystem::current_path();

    if (argc != 6)
    {
        std::cerr
            << "Usage: sbs-viewer.exe <scene specification json file> "
               "<path/to/vertex_shader.vs> <path/to/fragment_shader.fs> "
               "<path/to/wireframe_vertex_shader.vs> <path/to/wireframe_fragment_shader.fs>\n";
        return 1;
    }

    std::filesystem::path const scene_specification_path{argv[1]};
    std::filesystem::path const vertex_shader_path{argv[2]};
    std::filesystem::path const fragment_shader_path{argv[3]};
    std::filesystem::path const wireframe_vertex_shader_path{argv[4]};
    std::filesystem::path const wireframe_fragment_shader_path{argv[5]};

    sbs::rendering::renderer_t renderer{};
    sbs::physics::xpbd::solver_t solver{};

    std::vector<std::shared_ptr<sbs::common::shared_vertex_surface_mesh_i>>
        physical_node_surfaces{};

    double force_magnitude = 1000.;

    renderer.on_scene_loaded = [&](sbs::common::scene_t& scene) {
        solver.setup(scene.nodes);
        solver.timestep()   = 0.0167;
        solver.substeps()   = 30u;
        solver.iterations() = 30u;

        physical_node_surfaces.clear();
        for (auto const& node : scene.nodes)
        {
            if (!node->is_physically_simulated_body())
                continue;

            auto const tet_mesh =
                std::dynamic_pointer_cast<sbs::physics::xpbd::tetrahedral_mesh_t>(node);

            auto const phys_node =
                std::make_shared<sbs::physics::tetrahedral_mesh_surface_mesh_adapter_t>(
                    tet_mesh.get());
            physical_node_surfaces.push_back(phys_node);
        }

        renderer.pickers.clear();

        sbs::rendering::picker_t force_picker{&renderer, physical_node_surfaces};
        force_picker.should_picking_start = [](int button, int action, int mods) {
            return (
                button == GLFW_MOUSE_BUTTON_LEFT && mods == GLFW_MOD_SHIFT && action == GLFW_PRESS);
        };
        force_picker.should_picking_stop = [](bool left_mouse_button_pressed,
                                              bool right_mouse_button_pressed,
                                              bool middle_mouse_button_pressed,
                                              bool alt_key_pressed,
                                              bool ctrl_key_pressed,
                                              bool shift_key_pressed) {
            return !shift_key_pressed;
        };
        force_picker.should_pick =
            [](std::shared_ptr<sbs::common::shared_vertex_surface_mesh_i> node) {
                return true;
            };
        force_picker
            .on_picker_moved = [&force_magnitude](
                                   Eigen::Vector3d const& d,
                                   std::shared_ptr<sbs::common::shared_vertex_surface_mesh_i> node,
                                   std::uint32_t vi) {
            auto const phys_node =
                std::dynamic_pointer_cast<sbs::physics::tetrahedral_mesh_surface_mesh_adapter_t>(
                    node);
            auto const& tet_mesh_vertex_indices =
                phys_node->surface_to_tetrahedral_mesh_index_map();
            auto const tvi      = tet_mesh_vertex_indices[vi];
            auto const tet_mesh = phys_node->tetrahedral_mesh();
            tet_mesh->vertices().at(tvi).force() += d.normalized() * force_magnitude;
        };

        sbs::rendering::picker_t fix_picker{&renderer, physical_node_surfaces};
        fix_picker.should_picking_start = [](int button, int action, int mods) {
            return (
                button == GLFW_MOUSE_BUTTON_LEFT && mods == GLFW_MOD_CONTROL &&
                action == GLFW_PRESS);
        };
        fix_picker.should_picking_stop = [](bool left_mouse_button_pressed,
                                            bool right_mouse_button_pressed,
                                            bool middle_mouse_button_pressed,
                                            bool alt_key_pressed,
                                            bool ctrl_key_pressed,
                                            bool shift_key_pressed) {
            return !ctrl_key_pressed;
        };
        fix_picker.should_pick =
            [](std::shared_ptr<sbs::common::shared_vertex_surface_mesh_i> node) {
                return true;
            };
        fix_picker.picked = [](std::shared_ptr<sbs::common::shared_vertex_surface_mesh_i> node,
                               std::uint32_t vi) {
            auto const phys_node =
                std::dynamic_pointer_cast<sbs::physics::tetrahedral_mesh_surface_mesh_adapter_t>(
                    node);
            auto const& tet_mesh_vertex_indices =
                phys_node->surface_to_tetrahedral_mesh_index_map();
            auto const tvi                       = tet_mesh_vertex_indices[vi];
            auto const tet_mesh                  = phys_node->tetrahedral_mesh();
            tet_mesh->vertices().at(tvi).fixed() = !tet_mesh->vertices().at(tvi).fixed();
        };

        renderer.pickers.push_back(force_picker);
        renderer.pickers.push_back(fix_picker);
    };

    auto const step = [&](sbs::rendering::physics_timestep_throttler_t& throttler,
                          double frame_dt,
                          sbs::common::scene_t& scene) {
        throttler.timestep() = solver.timestep();
        solver.step();
        for (auto const& body : solver.simulated_bodies())
            body->mark_vertices_dirty();
    };
    sbs::rendering::physics_timestep_throttler_t throttler{solver.timestep(), step};
    renderer.on_new_physics_timestep = [&](double frame_dt, sbs::common::scene_t& scene) {
        throttler(frame_dt, scene);
    };

    auto const environment_node_factory = [](sbs::io::scene::scene_body_info const& sbi)
        -> std::shared_ptr<sbs::common::renderable_node_t> {
        return std::make_shared<sbs::common::static_mesh>(sbi.geometry);
    };

    auto const physics_node_factory = [](sbs::io::scene::physics_body_info const& pbi)
        -> std::shared_ptr<sbs::common::renderable_node_t> {
        sbs::physics::build_topology_parameters_t topological_params{};
        topological_params.edge_to_tetrahedra      = true;
        topological_params.tetrahedron_to_edge     = true;
        topological_params.triangle_to_tetrahedra  = true;
        topological_params.tetrahedron_to_triangle = true;

        sbs::physics::xpbd::simulation_parameters_t simulation_params{};

        auto node = std::make_shared<sbs::physics::xpbd::tetrahedral_mesh_t>(
            pbi.geometry,
            simulation_params,
            topological_params);

        for (sbs::physics::vertex_t& vertex : node->vertices())
        {
            vertex.velocity() = sbs::physics::vertex_t::velocity_type{
                pbi.velocity.vx,
                pbi.velocity.vy,
                pbi.velocity.vz};

            vertex.mass() = pbi.mass_density;
        }

        for (sbs::physics::tetrahedron_t& tetrahedron : node->tetrahedra())
        {
            tetrahedron.mass_density() = pbi.mass_density;
        }

        return node;
    };

    auto scene = sbs::io::load_scene(
        scene_specification_path,
        environment_node_factory,
        physics_node_factory);

    /*auto const cutting_surface_node_it = std::find_if(
        scene.nodes.begin(),
        scene.nodes.end(),
        [](std::shared_ptr<sbs::common::renderable_node_t> const object) {
            return object->id() == "cutter";
        });

    std::unique_ptr<sbs::physics::cutting::virtual_scalpel_h> virtual_scalpel{};
    if (cutting_surface_node_it != scene.nodes.end())
    {
        std::shared_ptr<sbs::common::renderable_node_t> const cutting_surface_node =
            *cutting_surface_node_it;

        auto const& p1 = cutting_surface_node->get_position(0u);
        auto const& p2 = cutting_surface_node->get_position(1u);

        virtual_scalpel = std::make_unique<sbs::physics::cutting::virtual_scalpel_h>(
            sbs::common::line_segment_t{
                sbs::common::point_t{p1.x, p1.y, p1.z},
                sbs::common::point_t{p2.x, p2.y, p2.z}}, );
    }*/

    std::shared_ptr<sbs::common::renderable_node_t> selected_node{};
    renderer.on_mouse_button_pressed = [&](GLFWwindow* window, int button, int action, int mods) {

    };

    // double xprev = 0., yprev = 0., dx = 0., dy = 0.;
    // bool is_first_mouse_move_callback           = true;
    // bool should_handle_cutting_surface_rotation = false;
    renderer.on_mouse_moved = [&](GLFWwindow* window, double x, double y) -> bool {
        // bool result = false;

        // int left_ctrl_key_state     = glfwGetKey(window, GLFW_KEY_LEFT_CONTROL);
        // int left_mouse_button_state = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);

        // if (left_ctrl_key_state == GLFW_PRESS && left_mouse_button_state == GLFW_PRESS)
        // {
        //     dx                                     = x - xprev;
        //     dy                                     = yprev - y;
        //     should_handle_cutting_surface_rotation = true;
        //     result                                 = true;
        // }

        // xprev = x;
        // yprev = y;
        return false;
    };

    renderer.on_new_imgui_frame = [&](sbs::common::scene_t& scene) {
        ImGui::Begin("Soft Body Simulator");

        static int selected_idx = 0;
        selected_node           = scene.nodes[selected_idx];
        if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_None))
        {
            ImGui::TreePush();
            float const scene_panel_width = ImGui::GetColumnWidth();
            if (ImGui::CollapsingHeader("Select Nodes##Scene", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ImGui::BulletText("Select active object");
                for (std::size_t b = 0u; b < scene.nodes.size(); ++b)
                {
                    auto const& body = scene.nodes[b];
                    ImGui::RadioButton(body->id().c_str(), &selected_idx, static_cast<int>(b));
                }

                selected_node = scene.nodes[selected_idx];
            }
            ImGui::TreePop();

            if (ImGui::Button("Reload", ImVec2(scene_panel_width / 2.f, 0.f)))
            {
                renderer.unload_current_scene();
                scene = sbs::io::load_scene(
                    scene_specification_path,
                    environment_node_factory,
                    physics_node_factory);
                renderer.load_scene(scene);
                selected_idx  = 0;
                selected_node = scene.nodes[selected_idx];
            }
            static bool should_render_wireframe = false;
            ImGui::Checkbox("Wireframe", &should_render_wireframe);
            if (should_render_wireframe)
            {
                selected_node->mark_should_render_wireframe();
            }
            else
            {
                selected_node->mark_should_render_triangles();
            }
        }

        if (ImGui::CollapsingHeader("Physics", ImGuiTreeNodeFlags_DefaultOpen))
        {
            bool const is_selected_node_a_physics_node =
                selected_node->is_physically_simulated_body();

            auto const tet_mesh =
                std::dynamic_pointer_cast<sbs::physics::xpbd::tetrahedral_mesh_t>(selected_node);

            ImGui::TreePush();
            ImGui::InputDouble("Fext magnitude", &force_magnitude, 100., 1000., "%.0f");
            ImGui::InputDouble(
                "Collision compliance",
                &solver.collision_compliance(),
                0.000001,
                0.001,
                "%.10f");

            if (ImGui::CollapsingHeader("XPBD") && is_selected_node_a_physics_node)
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
                static float young_modulus = 1e6f;
                static float poisson_ratio = 0.45f;
                ImGui::TreePush();
                ImGui::InputFloat("Young modulus##XPBD", &young_modulus, 1000.f, 10'000.f, "%.1f");
                ImGui::InputFloat("Poisson ratio##XPBD", &poisson_ratio, 0.01f, 0.1f, "%.2f");
                ImGui::TreePop();

                ImGui::Text("Solver");
                ImGui::TreePush();
                static float timestep = static_cast<float>(solver.timestep());
                static int iterations = static_cast<int>(solver.iterations());
                static int substeps   = static_cast<int>(solver.substeps());
                ImGui::InputFloat("Timestep##XPBD", &timestep, 0.001f, 0.01f);
                ImGui::InputInt("Iterations##XPBD", &iterations);
                ImGui::InputInt("Substeps##XPBD", &substeps);
                solver.timestep()   = static_cast<double>(timestep);
                solver.iterations() = static_cast<std::uint32_t>(iterations);
                solver.substeps()   = static_cast<std::uint32_t>(substeps);
                ImGui::TreePop();

                float const w = ImGui::GetColumnWidth();
                if (ImGui::Button("Apply##XPBD", ImVec2(w / 2.f, 0.f)))
                {
                    tet_mesh->simulation_parameters().alpha = static_cast<double>(alpha);
                    tet_mesh->simulation_parameters().hooke_coefficient =
                        1. / static_cast<double>(alpha);
                    tet_mesh->simulation_parameters().constraint_type =
                        constraint_type_choice == 0 ?
                            sbs::physics::xpbd::constraint_type_t::green :
                            sbs::physics::xpbd::constraint_type_t::distance;
                    tet_mesh->simulation_parameters().young_modulus =
                        static_cast<double>(young_modulus);
                    tet_mesh->simulation_parameters().poisson_ratio =
                        static_cast<double>(poisson_ratio);

                    renderer.unload_current_scene();
                    renderer.load_scene(scene);
                }
            }
            if (ImGui::CollapsingHeader("Mesh") && is_selected_node_a_physics_node)
            {
                std::size_t const num_vertices = tet_mesh->vertices().size();
                ImGui::Text("Vertices: %d", num_vertices);
                std::size_t const num_edges = tet_mesh->edges().size();
                ImGui::Text("Edges: %d", num_edges);
                std::size_t const num_triangles = tet_mesh->triangles().size();
                ImGui::Text("Triangles: %d", num_triangles);
                std::size_t const num_elements = tet_mesh->tetrahedra().size();
                ImGui::Text("Tetrahedra: %d", num_elements);
            }
            ImGui::TreePop();
            static bool are_physics_active = throttler.are_physics_active();
            ImGui::Checkbox("Activate physics", &are_physics_active);
            if (are_physics_active)
            {
                throttler.activate_physics();
            }
            else
            {
                throttler.deactivate_physics();
            }
            std::string const fps_str = "FPS: " + std::to_string(throttler.fps());
            ImGui::Text(fps_str.c_str());
        }

        /*if (ImGui::CollapsingHeader("Cutting", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (active_environment_node->id == "cutter")
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

                virtual_scalpel->set_translation(translation);
                if (should_handle_cutting_surface_rotation)
                {
                    virtual_scalpel->set_rotation_speed(static_cast<double>(sensitivity));
                    virtual_scalpel->set_yaw_axis(Eigen::Vector3d{
                        renderer.camera().front().x,
                        renderer.camera().front().y,
                        renderer.camera().front().z});
                    virtual_scalpel->set_pitch_axis(Eigen::Vector3d{
                        -renderer.camera().right().x,
                        -renderer.camera().right().y,
                        -renderer.camera().right().z});
                    virtual_scalpel->rotate(dx, dy);

                    should_handle_cutting_surface_rotation = false;
                }

                active_environment_node->render_model = virtual_scalpel->get_render_swept_surface();
                active_environment_node->render_state.should_transfer_vertices = true;

                ImGui::TreePop();
            }
        }*/

        ImGui::End();
    };

    bool const initialization_success = renderer.initialize();
    bool const shader_loading_success =
        renderer.use_shaders(vertex_shader_path, fragment_shader_path) &&
        renderer.use_wireframe_shaders(
            wireframe_vertex_shader_path,
            wireframe_fragment_shader_path);

    if (initialization_success && shader_loading_success)
    {
        renderer.load_scene(scene);
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