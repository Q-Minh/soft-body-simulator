#include <imgui/imgui.h>
#include <implot/implot.h>
#include <sbs/geometry/get_simple_bar_model.h>
#include <sbs/geometry/get_simple_plane_model.h>
#include <sbs/physics/body/environment_body.h>
#include <sbs/physics/collision/brute_force_cd_system.h>
#include <sbs/physics/mechanics/meshless_sph_body.h>
#include <sbs/physics/timestep.h>
#include <sbs/physics/xpbd/gauss_seidel_solver.h>
#include <sbs/physics/xpbd/meshless_sph_positional_constraint.h>
#include <sbs/physics/xpbd/meshless_sph_stvk_constraint.h>
#include <sbs/physics/xpbd/meshless_sph_surface_contact_handler.h>
#include <sbs/physics/xpbd/simulation.h>
#include <sbs/rendering/physics_timestep_throttler.h>
#include <sbs/rendering/pick.h>
#include <sbs/rendering/renderer.h>

int main(int argc, char** argv)
{
    /**
     * Setup simulation
     */
    sbs::physics::xpbd::simulation_t simulation{};
    simulation.simulation_parameters().compliance                  = 1e-10;
    simulation.simulation_parameters().damping                     = 1e-2;
    simulation.simulation_parameters().collision_compliance        = 1e-4;
    simulation.simulation_parameters().collision_damping           = 1e-5;
    simulation.simulation_parameters().poisson_ratio               = 0.45;
    simulation.simulation_parameters().young_modulus               = 1e6;
    simulation.simulation_parameters().positional_penalty_strength = 4.;

    sbs::common::geometry_t beam_geometry = sbs::geometry::get_simple_bar_model(12u, 4u, 12u);
    beam_geometry.set_color(255, 255, 0);
    sbs::scalar_type constexpr h = 2.;
    auto const beam_idx          = static_cast<sbs::index_type>(simulation.bodies().size());
    simulation.add_body();
    std::array<unsigned int, 3u> const particle_grid_resolution{12u, 4u, 12u};
    simulation.bodies()[beam_idx] = std::make_unique<sbs::physics::mechanics::meshless_sph_body_t>(
        simulation,
        beam_idx,
        beam_geometry,
        h,
        particle_grid_resolution);

    sbs::physics::mechanics::meshless_sph_body_t& beam =
        *dynamic_cast<sbs::physics::mechanics::meshless_sph_body_t*>(
            simulation.bodies()[beam_idx].get());
    Eigen::Affine3d beam_transform{Eigen::Translation3d(-3., 4., 2.)};
    beam_transform.rotate(
        Eigen::AngleAxisd(3.14159 / 2., Eigen::Vector3d{0., 1., 0.2}.normalized()));
    beam_transform.scale(Eigen::Vector3d{1, 0.4, 1});
    beam.transform(beam_transform);
    beam.initialize_physical_model();
    beam.initialize_visual_model();
    beam.initialize_collision_model();

    auto constexpr mass_density = 1.;
    sbs::scalar_type total_mass{0.};
    for (auto const& meshless_node : beam.nodes())
    {
        sbs::physics::xpbd::particle_t p{meshless_node.Xi()};
        p.mass() = mass_density * meshless_node.Vi();
        total_mass += p.mass();
        simulation.add_particle(p, beam_idx);
    }
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
    simulation.add_body(std::make_unique<sbs::physics::body::environment_body_t>(
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
        [](std::unique_ptr<sbs::physics::body::body_t>& b) { return &(b->collision_model()); });
    simulation.use_collision_detection_system(
        std::make_unique<sbs::physics::collision::brute_force_cd_system_t>(collision_objects));
    simulation.collision_detection_system()->use_contact_handler(
        std::make_unique<sbs::physics::xpbd::meshless_sph_surface_contact_handler_t>(simulation));

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
    timestep.solver()     = std::make_unique<sbs::physics::xpbd::gauss_seidel_solver_t>();

    sbs::rendering::physics_timestep_throttler_t throttler(
        timestep.dt(),
        [&](sbs::rendering::physics_timestep_throttler_t& throttler,
            double dt,
            sbs::physics::xpbd::simulation_t& simulation) { timestep.step(simulation); });

    renderer.on_new_physics_timestep = throttler;
    renderer.camera().position().x   = 0.;
    renderer.camera().position().y   = 5.;
    renderer.camera().position().z   = 25.;

    renderer.pickers.clear();
    std::vector<sbs::common::shared_vertex_surface_mesh_i*> surfaces_to_pick{};
    surfaces_to_pick.push_back(
        reinterpret_cast<sbs::common::shared_vertex_surface_mesh_i*>(&beam.surface_mesh()));

    std::unordered_map<sbs::index_type, Eigen::Vector3d> fixed_vertices{};
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
        if (fixed_vertices.find(vi) != fixed_vertices.end())
            return;

        auto* surface_mesh =
            reinterpret_cast<sbs::physics::mechanics::meshless_sph_surface_t*>(node);
        auto* mechanical_model = surface_mesh->mechanical_model();

        sbs::physics::mechanics::meshless_sph_surface_vertex_t const& meshless_surface_vertex =
            surface_mesh->embedded_surface_vertices()[vi];
        Eigen::Vector3d const target_position   = surface_mesh->vertex(vi).position;
        Eigen::Vector3d const& current_position = target_position;

        auto positional_constraint =
            std::make_unique<sbs::physics::xpbd::meshless_sph_positional_constraint_t>(
                simulation.simulation_parameters().collision_compliance,
                simulation.simulation_parameters().collision_damping,
                beam_idx,
                mechanical_model,
                meshless_surface_vertex,
                current_position,
                simulation.simulation_parameters().positional_penalty_strength,
                target_position);
        simulation.add_constraint(std::move(positional_constraint));

        fixed_vertices.insert({vi, target_position});
    };

    renderer.pickers.push_back(fix_picker);
    renderer.on_new_imgui_frame = [&](sbs::physics::xpbd::simulation_t& s) {
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
            std::string const element_count =
                "Elements: " + std::to_string(beam.topology().tetrahedron_count());
            ImGui::Text(element_count.c_str());
            std::string const total_mass_str = "Total mass: " + std::to_string(total_mass) + " g";
            ImGui::Text(total_mass_str.c_str());
        }

        static std::vector<sbs::scalar_type> Eis{};
        Eis.clear();
        Eis.reserve(beam.nodes().size());
        if (ImPlot::BeginPlot("Strains"))
        {
            for (auto const& meshless_node : beam.nodes())
            {
                Eigen::Matrix3d const& Ei    = meshless_node.Ei();
                sbs::scalar_type const Fnorm = Ei.squaredNorm();
                Eis.push_back(Fnorm);
            }
            ImPlot::PlotBars("||Ei||^2", Eis.data(), Eis.size());
            ImPlot::EndPlot();
        }

        ImGui::End();
    };

    renderer.on_pre_render = [&](sbs::physics::xpbd::simulation_t& s) {
        renderer.clear_points();
        for (auto const [vi, p] : fixed_vertices)
        {
            std::array<float, 9u> const vertex_attributes{
                static_cast<float>(p.x()),
                static_cast<float>(p.y()),
                static_cast<float>(p.z()),
                0.f,
                0.f,
                0.f,
                1.f,
                0.f,
                0.f};
            renderer.add_point(vertex_attributes);
        }
    };

    renderer.launch();

    return 0;
}