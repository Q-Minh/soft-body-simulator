#include <algorithm>
#include <array>
#include <imgui/imgui.h>
#include <implot/implot.h>
#include <iostream>
#include <sbs/geometry/get_simple_bar_model.h>
#include <sbs/geometry/get_simple_plane_model.h>
#include <sbs/math/mapping.h>
#include <sbs/physics/body/environment_body.h>
#include <sbs/physics/body/meshless_body.h>
#include <sbs/physics/collision/brute_force_cd_system.h>
#include <sbs/physics/mechanics/sph_meshless_model.h>
#include <sbs/physics/xpbd/contact_handler.h>
#include <sbs/physics/xpbd/gauss_seidel_solver.h>
#include <sbs/physics/xpbd/simulation.h>
#include <sbs/physics/xpbd/sph.h>
#include <sbs/physics/xpbd/timestep.h>
#include <sbs/rendering/physics_timestep_throttler.h>
#include <sbs/rendering/pick.h>
#include <sbs/rendering/renderer.h>

int main(int argc, char** argv)
{
    /**
     * Setup simulation
     */
    sbs::physics::xpbd::simulation_t simulation{};
    simulation.simulation_parameters().compliance                  = 1e-6;
    simulation.simulation_parameters().damping                     = 1e-4;
    simulation.simulation_parameters().collision_compliance        = 1e-4;
    simulation.simulation_parameters().collision_damping           = 1e-2;
    simulation.simulation_parameters().poisson_ratio               = 0.45;
    simulation.simulation_parameters().young_modulus               = 1e6;
    simulation.simulation_parameters().positional_penalty_strength = 4.;

    // Load geometry
    sbs::common::geometry_t beam_geometry = sbs::geometry::get_simple_bar_model(12u, 4u, 12u);
    beam_geometry.set_color(255, 255, 0);
    Eigen::Affine3d beam_transform{Eigen::Translation3d(-1., 4., 2.)};
    // beam_transform.rotate(
    //     Eigen::AngleAxisd(3.14159 / 2., Eigen::Vector3d{0., 1., 0.2}.normalized()));
    beam_transform.scale(Eigen::Vector3d{1., 0.4, 1.});
    beam_geometry                  = sbs::common::transform(beam_geometry, beam_transform);
    sbs::scalar_type const support = 1.;
    std::array<unsigned int, 3u> const resolution{12u, 4u, 12u};

    // Initialize soft body
    using kernel_function_type = sbs::math::poly6_kernel_t;
    using meshless_model_type = sbs::physics::mechanics::sph_meshless_model_t<kernel_function_type>;
    using visual_model_type =
        sbs::physics::visual::meshless_embedded_surface_t<meshless_model_type>;
    using body_type = sbs::physics::body::meshless_body_t<meshless_model_type>;

    auto const beam_idx = simulation.add_body();
    simulation.bodies()[beam_idx] =
        std::make_unique<body_type>(simulation, beam_idx, beam_geometry, resolution, support);
    body_type& beam = *dynamic_cast<body_type*>(simulation.bodies()[beam_idx].get());

    meshless_model_type& mechanical_model = beam.get_mechanical_model();
    sbs::scalar_type const mass_density   = 3.;
    for (auto i = 0u; i < mechanical_model.dof_count(); ++i)
    {
        Eigen::Vector3d const& Xi = mechanical_model.dof(i).cast<sbs::scalar_type>();
        sbs::physics::xpbd::particle_t p{Xi};
        sbs::index_type const ti   = mechanical_model.englobing_tetrahedron_of_particle(i);
        auto const N               = mechanical_model.particles_in_tetrahedron(ti).size();
        sbs::scalar_type const det = static_cast<sbs::scalar_type>(
            mechanical_model.domain().barycentric_map(ti).determinant());
        sbs::scalar_type const Vtet = (1. / 6.) * det;
        sbs::scalar_type const Vi   = Vtet / static_cast<sbs::scalar_type>(N);
        // sbs::scalar_type const Vi = mechanical_model.V(i);
        // p.mass() = mass_density * Vi;
        p.mass() = 1.;
        simulation.add_particle(p, beam_idx);
    }
    beam.get_visual_model().update();

    for (sbs::index_type i = 0u; i < mechanical_model.dof_count(); ++i)
    {
        auto const alpha = simulation.simulation_parameters().compliance;
        auto const beta  = simulation.simulation_parameters().damping;
        auto const nu    = simulation.simulation_parameters().poisson_ratio;
        auto const E     = simulation.simulation_parameters().young_modulus;

        // Use direct nodal integration with Shepard coefficient as "volume"
        // sbs::scalar_type const Vi   = mechanical_model.V(i);

        // Use direct nodal integration with uniform volume based on particle
        // sampling in tetrahedron. If N particles shared the same tet, then
        // the volume of each particle is Vi = Volume(tet) / N.
        sbs::index_type const ti   = mechanical_model.englobing_tetrahedron_of_particle(i);
        auto const N               = mechanical_model.particles_in_tetrahedron(ti).size();
        sbs::scalar_type const det = static_cast<sbs::scalar_type>(
            mechanical_model.domain().barycentric_map(ti).determinant());
        sbs::scalar_type const Vtet = (1. / 6.) * det;
        sbs::scalar_type const Vi   = Vtet / static_cast<sbs::scalar_type>(N);

        auto constraint =
            std::make_unique<sbs::physics::xpbd::stvk_sph_nodal_integration_constraint_t<
                meshless_model_type>>(alpha, beta, i, beam_idx, Vi, mechanical_model, E, nu);
        simulation.add_constraint(std::move(constraint));
    }

    beam.get_visual_model().set_color({1.f, 1.f, 0.f});

    // Initialize colliding rigid floor
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
    auto contact_handler = std::make_unique<sbs::physics::xpbd::contact_handler_t>(simulation);
    contact_handler->on_mesh_vertex_to_sdf_contact =
        [&](sbs::physics::collision::surface_mesh_particle_to_sdf_contact_t const& contact) {
            assert(contact.b1() == floor_idx);
            assert(contact.b2() == beam_idx);

            visual_model_type& surface = beam.get_visual_model();
            sbs::index_type const vi   = contact.vi();

            using interpolation_op_type = typename visual_model_type::interpolation_function_type;

            auto const alpha = simulation.simulation_parameters().collision_compliance;
            auto const beta  = simulation.simulation_parameters().collision_damping;

            using collision_constraint_type =
                sbs::physics::xpbd::sph_collision_constraint_t<visual_model_type>;
            auto collision_constraint = std::make_unique<collision_constraint_type>(
                alpha,
                beta,
                vi,
                beam_idx,
                surface,
                contact.point(),
                contact.normal());

            simulation.add_collision_constraint(std::move(collision_constraint));
        };
    simulation.collision_detection_system()->use_contact_handler(std::move(contact_handler));

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
    sbs::physics::xpbd::timestep_t timestep{};
    timestep.dt()         = 0.016;
    timestep.iterations() = 5u;
    timestep.substeps()   = 2u;
    timestep.solver()     = std::make_unique<sbs::physics::xpbd::gauss_seidel_solver_t>();

    sbs::rendering::physics_timestep_throttler_t throttler(
        timestep.dt(),
        [&](sbs::rendering::physics_timestep_throttler_t& throttler,
            double dt,
            sbs::physics::xpbd::simulation_t& simulation) { timestep.step(simulation); });

    renderer.on_new_physics_timestep = throttler;
    renderer.camera().position().x   = 0.;
    renderer.camera().position().y   = 5.;
    renderer.camera().position().z   = 40.;

    renderer.pickers.clear();
    /*std::vector<sbs::common::shared_vertex_surface_mesh_i*> surfaces_to_pick{};
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
        auto const tvi                    = tet_mesh_boundary->from_surface_vertex(vi);
        auto const tet_mesh               = tet_mesh_boundary->tetrahedral_mesh();
        sbs::physics::xpbd::particle_t& p = simulation.particles()[beam_idx][tvi];
        p.mass()                          = p.fixed() ? 1. : 0.;
    };

    renderer.pickers.push_back(fix_picker);*/
    bool should_render_points{false};
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
            ImGui::Checkbox(
                "Render points##Scene",
                [&]() { return should_render_points; },
                [&](bool value) { should_render_points = value; });
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

            if (ImGui::Button("Step##Physics"))
            {
                timestep.step(simulation);
            }
            ImGui::Button("Hold Step##Physics");
            if (ImGui::IsItemActive())
            {
                timestep.step(simulation);
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

            std::size_t const num_dofs     = mechanical_model.dof_count();
            std::string const num_dofs_str = "DOFs: " + std::to_string(num_dofs);
            ImGui::Text(num_dofs_str.c_str());

            std::size_t const num_constraints = simulation.constraints().size();
            std::string const num_constraints_str =
                "Constraints: " + std::to_string(num_constraints);
            ImGui::Text(num_constraints_str.c_str());
        }

        ImGui::End();
    };

    renderer.on_pre_render = [&](sbs::physics::xpbd::simulation_t& s) {
        if (should_render_points)
        {
            renderer.clear_points();
            auto const& particles = s.particles()[beam_idx];
            for (auto j = 0; j < mechanical_model.dof_count(); ++j)
            {
                auto const& p = particles[j];
                std::array<float, 3u> const vp{
                    static_cast<float>(p.x().x()),
                    static_cast<float>(p.x().y()),
                    static_cast<float>(p.x().z())};
                std::array<float, 3u> const vn{0.f, 0.f, 0.f};
                std::array<float, 3u> const vc{1.f, 0.f, 0.f};
                renderer.add_point(vp, vn, vc);
            }
        }
    };

    renderer.launch();

    return 0;
}