#include <imgui/imgui.h>
#include <sbs/geometry/get_simple_bar_model.h>
#include <sbs/geometry/get_simple_plane_model.h>
#include <sbs/physics/collision/brute_force_cd_system.h>
#include <sbs/physics/environment_body.h>
#include <sbs/physics/gauss_seidel_solver.h>
#include <sbs/physics/mechanics/meshless_sph_body.h>
#include <sbs/physics/simulation.h>
#include <sbs/physics/timestep.h>
#include <sbs/physics/xpbd/contact_handler.h>
#include <sbs/physics/xpbd/meshless_sph_stvk_constraint.h>
#include <sbs/rendering/physics_timestep_throttler.h>
#include <sbs/rendering/pick.h>
#include <sbs/rendering/renderer.h>

int main(int argc, char** argv)
{
    /**
     * Setup simulation
     */
    sbs::physics::simulation_t simulation{};
    simulation.simulation_parameters().damping           = 1.;
    simulation.simulation_parameters().collision_damping = 0.01;
    simulation.simulation_parameters().poisson_ratio     = 0.3;
    simulation.simulation_parameters().young_modulus     = 1e6;

    sbs::common::geometry_t beam_geometry = sbs::geometry::get_simple_bar_model(5u, 5u, 20u);
    beam_geometry.set_color(255, 255, 0);
    sbs::scalar_type constexpr h = 1.;
    auto const beam_idx          = static_cast<sbs::index_type>(simulation.bodies().size());
    simulation.add_body();
    simulation.bodies()[beam_idx] = std::make_unique<sbs::physics::mechanics::meshless_sph_body_t>(
        simulation,
        beam_idx,
        beam_geometry,
        h);

    sbs::physics::mechanics::meshless_sph_body_t& beam =
        *dynamic_cast<sbs::physics::mechanics::meshless_sph_body_t*>(
            simulation.bodies()[beam_idx].get());
    Eigen::Affine3d beam_transform{Eigen::Translation3d(-3., 5., -1.)};
    beam_transform.rotate(
        Eigen::AngleAxisd(3.14159 / 2., Eigen::Vector3d{0., 1., 0.2}.normalized()));
    beam_transform.scale(Eigen::Vector3d{0.2, 0.2, 0.5});
    beam.transform(beam_transform);
    beam.initialize_physical_model();
    beam.initialize_visual_model();
    beam.initialize_collision_model();

    for (std::size_t i = 0u; i < beam.nodes().size(); ++i)
    {
        auto const alpha = simulation.simulation_parameters().compliance;
        auto const beta  = simulation.simulation_parameters().damping;
        auto const E     = simulation.simulation_parameters().young_modulus;
        auto const nu    = simulation.simulation_parameters().poisson_ratio;
        sbs::physics::mechanics::meshless_sph_node_t& node = beam.nodes()[i];
        auto const ni                                      = static_cast<sbs::index_type>(i);
        auto constraint = std::make_unique<sbs::physics::xpbd::meshless_sph_stvk_constraint_t>(
            alpha,
            beta,
            simulation,
            beam_idx,
            ni,
            E,
            nu,
            node);
        simulation.add_constraint(std::move(constraint));
    }

    sbs::common::geometry_t floor_geometry =
        sbs::geometry::get_simple_plane_model({-20., -20.}, {20., 20.}, 0., 1e-2);
    floor_geometry.set_color(100, 100, 100);
    auto const floor_idx = static_cast<sbs::index_type>(simulation.bodies().size());
    Eigen::AlignedBox3d const floor_volume{
        Eigen::Vector3d{-20., -5., -20.},
        Eigen::Vector3d{20., 5., 20.}};
    auto const floor_collision_model = sbs::physics::collision::sdf_model_t::from_plane(
        Eigen::Hyperplane<sbs::scalar_type, 3>(
            Eigen::Vector3d{0., 1., 0.},
            Eigen::Vector3d{0., 0., 0.}),
        floor_volume);
    simulation.add_body(std::make_unique<sbs::physics::environment_body_t>(
        simulation,
        floor_idx,
        floor_geometry,
        floor_collision_model));

    /**
     * Setup collision detection
     */
    std::vector<sbs::physics::collision::collision_model_t*> collision_objects{};
    std::transform(
        simulation.bodies().begin(),
        simulation.bodies().end(),
        std::back_inserter(collision_objects),
        [](std::unique_ptr<sbs::physics::body_t>& b) { return &(b->collision_model()); });
    simulation.use_collision_detection_system(
        std::make_unique<sbs::physics::collision::brute_force_cd_system_t>(collision_objects));
    simulation.collision_detection_system()->use_contact_handler(
        std::make_unique<sbs::physics::xpbd::contact_handler_t>(simulation));

    /**
     * Setup renderer
     */
    if (argc != 2)
    {
        std::cerr << "Usage: tester.exe <path/to/shader/directory/>\n";
        return 1;
    }

    sbs::rendering::point_light_t point_light(0.f, 5.f, 0.f);
    sbs::rendering::directional_light_t directional_light(0.f, 0.f, -1.f);
    sbs::rendering::renderer_t renderer(simulation, point_light, directional_light);

    std::filesystem::path const cwd = std::filesystem::current_path();

    std::filesystem::path const shader_directory_path{argv[1]};
    std::filesystem::path const vertex_shader_path   = shader_directory_path / "mesh.vs";
    std::filesystem::path const fragment_shader_path = shader_directory_path / "mesh.fs";
    std::filesystem::path const wireframe_vertex_shader_path =
        shader_directory_path / "wireframe.vs";
    std::filesystem::path const wireframe_fragment_shader_path =
        shader_directory_path / "wireframe.fs";
    std::filesystem::path const point_vertex_shader_path   = shader_directory_path / "point.vs";
    std::filesystem::path const point_fragment_shader_path = shader_directory_path / "point.fs";

    renderer.use_mesh_shaders(vertex_shader_path, fragment_shader_path);
    renderer.use_wireframe_shaders(wireframe_vertex_shader_path, wireframe_fragment_shader_path);
    renderer.use_point_shaders(point_vertex_shader_path, point_fragment_shader_path);

    /**
     * Setup time integration technique
     */
    sbs::physics::timestep_t timestep{};
    timestep.dt()         = 0.016;
    timestep.iterations() = 5u;
    timestep.substeps()   = 1u;
    timestep.solver()     = std::make_unique<sbs::physics::gauss_seidel_solver_t>();

    sbs::rendering::physics_timestep_throttler_t throttler(
        timestep.dt(),
        [&](sbs::rendering::physics_timestep_throttler_t& throttler,
            double dt,
            sbs::physics::simulation_t& simulation) { timestep.step(simulation); });

    renderer.on_new_physics_timestep = throttler;
    renderer.camera().position().x   = 0.;
    renderer.camera().position().y   = 5.;
    renderer.camera().position().z   = 25.;

    renderer.pickers.clear();
    std::vector<sbs::common::shared_vertex_surface_mesh_i*> surfaces_to_pick{};
    surfaces_to_pick.push_back(
        reinterpret_cast<sbs::common::shared_vertex_surface_mesh_i*>(&beam.surface_mesh()));

    sbs::rendering::picker_t fix_picker{&renderer, surfaces_to_pick};
    fix_picker.should_picking_start = [](int button, int action, int mods) {
        return (
            button == GLFW_MOUSE_BUTTON_LEFT && mods == GLFW_MOD_CONTROL && action == GLFW_PRESS);
    };
    fix_picker.should_picking_stop = [](bool left_mouse_button_pressed,
                                        bool right_mouse_button_pressed,
                                        bool middle_mouse_button_pressed,
                                        bool alt_key_pressed,
                                        bool ctrl_key_pressed,
                                        bool shift_key_pressed) {
        return !ctrl_key_pressed;
    };
    fix_picker.should_pick = [](sbs::common::shared_vertex_surface_mesh_i* node) {
        return true;
    };
    fix_picker.picked = [&](sbs::common::shared_vertex_surface_mesh_i* node, std::uint32_t vi) {
        auto* tet_mesh_boundary =
            reinterpret_cast<sbs::physics::tetrahedral_mesh_boundary_t*>(node);
        auto const tvi              = tet_mesh_boundary->from_surface_vertex(vi);
        auto const tet_mesh         = tet_mesh_boundary->tetrahedral_mesh();
        sbs::physics::particle_t& p = simulation.particles()[beam_idx][tvi];
        p.mass()                    = p.fixed() ? 1. : 0.;
    };

    renderer.pickers.push_back(fix_picker);
    renderer.on_new_imgui_frame = [&](sbs::physics::simulation_t& s) {
        ImGui::Begin("Soft Body Simulator");

        static int selected_idx = 0;
        if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_None))
        {
            ImGui::TreePush();
            float const scene_panel_width = ImGui::GetColumnWidth();
            if (ImGui::CollapsingHeader("Select Nodes##Scene", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ImGui::BulletText("Select active object");
                for (std::size_t b = 0u; b < s.bodies().size(); ++b)
                {
                    auto const& body = s.bodies()[b];
                    ImGui::RadioButton(
                        std::to_string(body->id()).c_str(),
                        &selected_idx,
                        static_cast<int>(b));
                }
            }
            ImGui::TreePop();
            static bool should_render_wireframe = false;
            ImGui::Checkbox("Wireframe", &should_render_wireframe);
            if (should_render_wireframe)
            {
                s.bodies()[selected_idx]->visual_model().mark_should_render_wireframe();
            }
            else
            {
                s.bodies()[selected_idx]->visual_model().mark_should_render_triangles();
            }
        }

        if (ImGui::CollapsingHeader("Physics", ImGuiTreeNodeFlags_DefaultOpen))
        {
            auto* throttler_ptr = renderer.on_new_physics_timestep
                                      .target<sbs::rendering::physics_timestep_throttler_t>();
            static bool are_physics_active = throttler_ptr->are_physics_active();
            ImGui::Checkbox("Activate physics", &are_physics_active);
            if (are_physics_active)
            {
                throttler_ptr->activate_physics();
            }
            else
            {
                throttler_ptr->deactivate_physics();
            }
            static std::size_t windowed_fps_sum     = 0u;
            static std::size_t num_frames           = 0u;
            static std::size_t windowed_average_fps = 0u;
            auto const fps                          = throttler_ptr->fps();
            windowed_fps_sum += fps;
            ++num_frames;

            std::size_t constexpr window_size = 50u;
            if (num_frames == window_size)
            {
                windowed_average_fps = windowed_fps_sum / num_frames;
                num_frames           = 0u;
                windowed_fps_sum     = 0u;
            }

            std::string const fps_str = "FPS: " + std::to_string(fps);
            ImGui::Text(fps_str.c_str());
            std::string const fps_avg_str = "Mean FPS: " + std::to_string(windowed_average_fps);
            ImGui::Text(fps_avg_str.c_str());
            std::string const particle_count = "Particles: " + std::to_string(beam.nodes().size());
            ImGui::Text(particle_count.c_str());
        }

        ImGui::End();
    };

    renderer.on_pre_render = [&](sbs::physics::simulation_t& s) {
        renderer.clear_points();
        for (auto const& particles : s.particles())
        {
            for (auto const& p : particles)
            {
                if (p.fixed())
                {
                    std::array<float, 9u> const vertex_attributes{
                        static_cast<float>(p.x().x()),
                        static_cast<float>(p.x().y()),
                        static_cast<float>(p.x().z()),
                        0.f,
                        0.f,
                        0.f,
                        1.f,
                        0.f,
                        0.f};
                    renderer.add_point(vertex_attributes);
                }
            }
        }
    };

    renderer.launch();

    return 0;
}