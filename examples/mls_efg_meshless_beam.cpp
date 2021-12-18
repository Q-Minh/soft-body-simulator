#include <algorithm>
#include <array>
#include <imgui/imgui.h>
#include <iostream>
#include <sbs/geometry/get_simple_bar_model.h>
#include <sbs/geometry/get_simple_plane_model.h>
#include <sbs/math/mapping.h>
#include <sbs/math/mls.h>
#include <sbs/math/quadrature.h>
#include <sbs/physics/body/environment_body.h>
#include <sbs/physics/body/meshless_body.h>
#include <sbs/physics/collision/brute_force_cd_system.h>
#include <sbs/physics/mechanics/efg_tetrahedral_meshless_model.h>
#include <sbs/physics/xpbd/contact_handler.h>
#include <sbs/physics/xpbd/gauss_seidel_solver.h>
#include <sbs/physics/xpbd/interpolated_particle_collision_constraint.h>
#include <sbs/physics/xpbd/simulation.h>
#include <sbs/physics/xpbd/strain_energy_constraint.h>
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

    // Load geometry
    sbs::common::geometry_t beam_geometry = sbs::geometry::get_simple_bar_model(4u, 4u, 4u);
    beam_geometry.set_color(255, 255, 0);
    //Eigen::Affine3d beam_transform{Eigen::Translation3d(-10., 5., -1.)};
    //beam_transform.rotate(
    //    Eigen::AngleAxisd(3.14159 / 2., Eigen::Vector3d{0., 1., 0.2}.normalized()));
    //beam_transform.scale(Eigen::Vector3d{1.0, 0.8, 2.});
    //beam_geometry                  = sbs::common::transform(beam_geometry, beam_transform);
    sbs::scalar_type const support = 1.3;
    std::array<unsigned int, 3u> const resolution{6u, 6u, 6u};

    // Initialize soft body
    unsigned int constexpr order = 1u;
    using kernel_function_type   = sbs::math::poly6_kernel_t;
    using meshless_model_type =
        sbs::physics::mechanics::efg_tetrahedral_meshless_model_t<kernel_function_type, order>;
    using body_type = sbs::physics::body::meshless_body_t<meshless_model_type>;

    auto const beam_idx = simulation.add_body();
    simulation.bodies()[beam_idx] =
        std::make_unique<body_type>(simulation, beam_idx, beam_geometry, resolution, support);
    body_type& beam = *dynamic_cast<body_type*>(simulation.bodies()[beam_idx].get());

    meshless_model_type& mechanical_model = beam.get_mechanical_model();
    for (auto i = 0u; i < mechanical_model.dof_count(); ++i)
    {
        Eigen::Vector3d const& Xi = mechanical_model.point(i).cast<sbs::scalar_type>();
        sbs::physics::xpbd::particle_t p{Xi};
        p.mass() = 1.;
        simulation.add_particle(p, beam_idx);
    }
    sbs::geometry::tetrahedral_domain_t const& domain = mechanical_model.domain();
    auto const& topology                              = domain.topology();
    auto const tetrahedron_count                      = topology.tetrahedron_count();
    for (sbs::index_type ti = 0u; ti < tetrahedron_count; ++ti)
    {
        sbs::topology::tetrahedron_t const& t = topology.tetrahedron(ti);
        Eigen::Vector3d const& p1             = domain.position(t.v1());
        Eigen::Vector3d const& p2             = domain.position(t.v2());
        Eigen::Vector3d const& p3             = domain.position(t.v3());
        Eigen::Vector3d const& p4             = domain.position(t.v4());

        // Insert integration point at the barycenter of the integration domain's tetrahedra

        sbs::math::tetrahedron_affine_mapping_t const mapping(p1, p2, p3, p4);
        sbs::math::tetrahedron_1point_constant_quadrature_rule_t<decltype(mapping)> quadrature_rule{
            mapping};

        for (auto i = 0u; i < quadrature_rule.points.size(); ++i)
        {
            auto const Xi = quadrature_rule.points[i];
            auto const wi = quadrature_rule.weights[i];

            Eigen::Vector3d const integration_point = Xi.cast<sbs::scalar_type>();
            sbs::index_type const integration_point_idx =
                static_cast<sbs::index_type>(mechanical_model.integration_point_count());
            bool const was_integration_point_added =
                mechanical_model.add_integration_point(integration_point);

            if (!was_integration_point_added)
                continue;

            using interpolation_op_type = typename meshless_model_type::interpolation_function_type;
            interpolation_op_type const interpolation_field_around_integration_point =
                mechanical_model.interpolation_field_from_integration_point(integration_point_idx);

            using constraint_type =
                sbs::physics::xpbd::strain_energy_quadrature_constraint_t<interpolation_op_type>;

            std::vector<sbs::index_type> const particle_indices =
                mechanical_model.neighbours_of_integration_point(integration_point_idx);

            std::vector<sbs::index_type> particle_body_indices(particle_indices.size(), beam_idx);

            auto const alpha = simulation.simulation_parameters().compliance;
            auto const beta  = simulation.simulation_parameters().damping;
            auto const nu    = simulation.simulation_parameters().poisson_ratio;
            auto const E     = simulation.simulation_parameters().young_modulus;

            auto constraint = std::make_unique<constraint_type>(
                alpha,
                beta,
                particle_indices,
                particle_body_indices,
                interpolation_field_around_integration_point,
                wi,
                Xi,
                E,
                nu);

            simulation.add_constraint(std::move(constraint));
        }
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

            using visual_model_type =
                sbs::physics::visual::meshless_embedded_surface_t<meshless_model_type>;

            visual_model_type const& surface = beam.get_visual_model();
            sbs::index_type const vi         = contact.vi();
            Eigen::Vector3d const& Xi        = surface.reference_position(vi);

            using interpolation_op_type = typename visual_model_type::interpolation_function_type;
            interpolation_op_type const interpolation_op = surface.interpolation_operator(vi);

            auto const alpha = simulation.simulation_parameters().compliance;
            auto const beta  = simulation.simulation_parameters().damping;

            // TODO: Implement this, but neighbours should be precomputed
            std::vector<sbs::index_type> const js = surface.meshless_neighbours_of_vertex(vi);
            std::vector<sbs::index_type> bis(js.size(), beam_idx);

            using constraint_type =
                sbs::physics::xpbd::interpolated_particle_collision_constraint_t<
                    interpolation_op_type>;
            auto collision_constraint = std::make_unique<constraint_type>(
                alpha,
                beta,
                js,
                bis,
                interpolation_op,
                Xi,
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
        }

        ImGui::End();
    };

    renderer.on_pre_render = [&](sbs::physics::xpbd::simulation_t& s) {
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