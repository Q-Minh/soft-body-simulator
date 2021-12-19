#include <algorithm>
#include <array>
#include <imgui/imgui.h>
#include <implot/implot.h>
#include <iostream>
#include <sbs/geometry/curves.h>
#include <sbs/geometry/get_simple_bar_model.h>
#include <sbs/geometry/get_simple_plane_model.h>
#include <sbs/geometry/transform.h>
#include <sbs/math/kernels.h>
#include <sbs/math/mapping.h>
#include <sbs/math/quadrature.h>
#include <sbs/physics/body/environment_body.h>
#include <sbs/physics/body/fem_mixed_body.h>
#include <sbs/physics/collision/brute_force_cd_system.h>
#include <sbs/physics/mechanics/fem_sph_efg_model.h>
#include <sbs/physics/xpbd/contact_handler.h>
#include <sbs/physics/xpbd/fem.h>
#include <sbs/physics/xpbd/fem_mixed_constraint.h>
#include <sbs/physics/xpbd/gauss_seidel_solver.h>
#include <sbs/physics/xpbd/simulation.h>
#include <sbs/physics/xpbd/timestep.h>
#include <sbs/rendering/physics_timestep_throttler.h>
#include <sbs/rendering/pick.h>
#include <sbs/rendering/renderer.h>

using kernel_function_type = sbs::math::poly6_kernel_t;
using fem_mixed_model_type = sbs::physics::mechanics::fem_sph_efg_model_t<kernel_function_type>;
using fem_model_type       = typename fem_mixed_model_type::fem_model_type;
using meshless_model_type  = typename fem_mixed_model_type::meshless_model_type;
using visual_model_type = sbs::physics::visual::fem_mixed_embedded_surface_t<fem_mixed_model_type>;
using body_type         = sbs::physics::body::fem_mixed_body_t<fem_mixed_model_type>;

std::size_t num_boundary_tets(fem_mixed_model_type const& body)
{
    auto const& topology    = body.domain().topology();
    std::size_t const count = std::count_if(
        topology.tetrahedra().begin(),
        topology.tetrahedra().end(),
        [&](sbs::topology::tetrahedron_t const& t) { return topology.is_boundary_tetrahedron(t); });
    return count;
}

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: <>.exe <path/to/shader/directory/>\n";
        return 1;
    }

    /**
     * Setup simulation
     */
    sbs::physics::xpbd::simulation_t simulation{};
    simulation.simulation_parameters().compliance                  = 1e-10;
    simulation.simulation_parameters().damping                     = 0.;
    simulation.simulation_parameters().collision_compliance        = 1e-6;
    simulation.simulation_parameters().collision_damping           = 1e-2;
    simulation.simulation_parameters().poisson_ratio               = 0.45;
    simulation.simulation_parameters().young_modulus               = 1e6;
    simulation.simulation_parameters().positional_penalty_strength = 4.;

    // Load geometry
    sbs::common::geometry_t beam_geometry = sbs::geometry::get_simple_bar_model(16u, 4u, 4u);
    beam_geometry.set_color(255, 255, 0);
    Eigen::Affine3d beam_transform{Eigen::Translation3d(-1., 4., 2.)};
    //beam_transform.rotate(
    //    Eigen::AngleAxisd(3.14159 / 16., Eigen::Vector3d{0., 0., 1.}.normalized()));
    // beam_transform.scale(Eigen::Vector3d{1., 0.4, 1.});
    beam_geometry                  = sbs::common::transform(beam_geometry, beam_transform);
    sbs::scalar_type const support = 1.1;
    std::array<unsigned int, 3u> const resolution{8u, 4u, 4u};

    // Initialize soft body
    auto const beam_idx = simulation.add_body();
    simulation.bodies()[beam_idx] =
        std::make_unique<body_type>(simulation, beam_idx, beam_geometry, resolution, support);
    body_type& beam = *dynamic_cast<body_type*>(simulation.bodies()[beam_idx].get());

    fem_mixed_model_type& mechanical_model = beam.get_mechanical_model();
    sbs::scalar_type const mass_density    = 3.;

    // Variables for Dirichlet boundary condition imposition
    std::vector<sbs::index_type> left_bs{};
    std::vector<sbs::index_type> right_bs{};
    std::vector<sbs::index_type> left_fixed_particles{};
    std::vector<Eigen::Vector3d> left_dirichlet_positions{};
    std::vector<sbs::index_type> right_fixed_particles{};
    std::vector<Eigen::Vector3d> right_dirichlet_positions{};

    // Create FEM particles
    auto const fem_particle_index_offset = simulation.particles()[beam_idx].size();
    for (auto i = 0u; i < mechanical_model.dof_count(); ++i)
    {
        Eigen::Vector3d const& xi = mechanical_model.dof(i);
        sbs::physics::xpbd::particle_t p{xi};
        p.mass() = 1.;
        sbs::index_type const idx =
            static_cast<sbs::index_type>(simulation.particles()[beam_idx].size());
        if (p.x0().x() < 0.)
        {
            left_fixed_particles.push_back(idx);
            left_dirichlet_positions.push_back(p.x0());
        }
        // if (p.x0().x() > 5.)
        //{
        //     right_fixed_particles.push_back(idx);
        //     right_dirichlet_positions.push_back(p.x0());
        // }
        simulation.add_particle(p, beam_idx);
    }

    // Create SPH particles
    auto const meshless_particle_index_offset = simulation.particles()[beam_idx].size();
    auto& meshless_model                      = mechanical_model.meshless_model();
    for (auto j = 0u; j < meshless_model.dof_count(); ++j)
    {
        Eigen::Vector3d const& xj = meshless_model.dof(j);
        sbs::physics::xpbd::particle_t p{xj};
        p.mass() = 1.;
        sbs::index_type const idx =
            static_cast<sbs::index_type>(simulation.particles()[beam_idx].size());
        if (p.x0().x() < 0.)
        {
            left_fixed_particles.push_back(idx);
            left_dirichlet_positions.push_back(p.x0());
        }
        // if (p.x0().x() > 5.)
        //{
        //     right_fixed_particles.push_back(idx);
        //     right_dirichlet_positions.push_back(p.x0());
        // }
        simulation.add_particle(p, beam_idx);
    }
    beam.set_fem_particle_offset(fem_particle_index_offset);
    beam.set_meshless_particle_offset(meshless_particle_index_offset);

    beam.get_visual_model().update();

    // Apply initial Dirichlet BCs
    left_bs.resize(left_fixed_particles.size(), beam_idx);
    right_bs.resize(right_fixed_particles.size(), beam_idx);

    simulation.apply_dirichlet_boundary_conditions(
        left_bs,
        left_fixed_particles,
        left_dirichlet_positions);

    simulation.apply_dirichlet_boundary_conditions(
        right_bs,
        right_fixed_particles,
        right_dirichlet_positions);

    // Create constraints
    auto const& domain   = mechanical_model.domain();
    auto const& topology = domain.topology();

    auto const alpha = simulation.simulation_parameters().compliance;
    auto const beta  = simulation.simulation_parameters().damping;
    auto const nu    = simulation.simulation_parameters().poisson_ratio;
    auto const E     = simulation.simulation_parameters().young_modulus;

    bool is_fully_mixed{true};
    std::size_t num_fully_meshed_cells{0u};
    std::size_t num_meshless_constraints{0u};
    std::size_t num_boundary_unmixed_tets{0u};
    std::size_t num_mixed_tets{0u};
    // EFG constraints
    for (auto k = 0u; k < mechanical_model.efg_integration_point_count(); ++k)
    {
        Eigen::Vector3d const& Xk = mechanical_model.integration_point(k);
        sbs::index_type const e   = mechanical_model.tetrahedron_around_integration_point(k);
        auto const& element       = mechanical_model.element(e);
        sbs::scalar_type const det =
            static_cast<sbs::scalar_type>(domain.barycentric_map(e).determinant());
        sbs::scalar_type const Vtet = (1. / 6.) * det;

        using constraint_type =
            sbs::physics::xpbd::stvk_fem_mixed_efg_integration_constraint_t<fem_mixed_model_type>;

        auto constraint = std::make_unique<constraint_type>(
            alpha,
            beta,
            k,
            beam_idx,
            Vtet,
            e,
            mechanical_model,
            E,
            nu,
            fem_particle_index_offset,
            meshless_particle_index_offset);

        simulation.add_constraint(std::move(constraint));
    }
    // FEM constraints
    for (auto e = 0u; e < mechanical_model.element_count(); ++e)
    {
        bool const is_mixed_tet = mechanical_model.is_mixed_cell(e);
        if (is_mixed_tet)
            continue;

        auto const& cell = mechanical_model.cell(e);

        auto const v1 = cell.node(0u);
        auto const v2 = cell.node(1u);
        auto const v3 = cell.node(2u);
        auto const v4 = cell.node(3u);

        auto const& element = mechanical_model.element(e);
        auto const& mapping = element.mapping();
        sbs::math::tetrahedron_1point_constant_quadrature_rule_t<decltype(mapping)> quadrature_rule{
            mapping};

        auto const& Xg = quadrature_rule.points.front().cast<sbs::scalar_type>();
        auto const& wg = static_cast<sbs::scalar_type>(quadrature_rule.weights.front());

        using constraint_type = sbs::physics::xpbd::stvk_tetrahedral_quadrature_strain_constraint_t<
            typename fem_mixed_model_type::cell_type>;

        std::vector<sbs::index_type> is{v1, v2, v3, v4};
        std::vector<sbs::index_type> bs{beam_idx, beam_idx, beam_idx, beam_idx};

        auto const& interpolation_function = mechanical_model.fem_interpolation_field_at(e);

        auto constraint = std::make_unique<constraint_type>(
            alpha,
            beta,
            fem_particle_index_offset,
            is,
            bs,
            interpolation_function,
            wg,
            Xg,
            E,
            nu);

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
                sbs::physics::xpbd::fem_mixed_collision_constraint_t<visual_model_type>;
            auto collision_constraint = std::make_unique<collision_constraint_type>(
                alpha,
                beta,
                vi,
                beam_idx,
                surface,
                contact.point(),
                contact.normal(),
                fem_particle_index_offset,
                meshless_particle_index_offset);

            simulation.add_collision_constraint(std::move(collision_constraint));
        };
    simulation.collision_detection_system()->use_contact_handler(std::move(contact_handler));

    /**
     * Setup renderer
     */
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

    std::vector<Eigen::Vector3d> control_points{};
    control_points.push_back(Eigen::Vector3d{0., 0., 0.});
    control_points.push_back(Eigen::Vector3d{1., 1., 0.});
    control_points.push_back(Eigen::Vector3d{2., 0., 0.});
    control_points.push_back(Eigen::Vector3d{2., 0., 0.});
    control_points.push_back(Eigen::Vector3d{3., -1., 0.});
    control_points.push_back(Eigen::Vector3d{4., 0., 0.});

    bool should_render_points{false};
    bool should_render_active_fem_nodes{false};
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
            ImGui::Checkbox(
                "Render active nodes##Scene",
                [&]() { return should_render_active_fem_nodes; },
                [&](bool value) { should_render_active_fem_nodes = value; });
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

            std::string const fem_dofs_str =
                "FEM dofs: " + std::to_string(mechanical_model.num_active_nodes());
            std::string const sph_dofs_str =
                "SPH dofs: " + std::to_string(mechanical_model.num_sph_nodes());
            ImGui::Text(fem_dofs_str.c_str());
            ImGui::Text(sph_dofs_str.c_str());

            std::size_t const num_tetrahedra =
                mechanical_model.domain().topology().tetrahedron_count();
            std::string const num_tetrahedra_str = "Tets: " + std::to_string(num_tetrahedra);
            ImGui::Text(num_tetrahedra_str.c_str());

            std::size_t const num_boundary_tetrahedra = num_boundary_tets(mechanical_model);
            std::string const num_boundary_tetrahedra_str =
                "Boundary tets: " + std::to_string(num_boundary_tetrahedra);
            ImGui::Text(num_boundary_tetrahedra_str.c_str());

            std::size_t const num_constraints = simulation.constraints().size();
            std::string const num_constraints_str =
                "Constraints: " + std::to_string(num_constraints);
            ImGui::Text(num_constraints_str.c_str());

            static sbs::scalar_type translation = 0.f;
            static float translation_speed      = 0.01f;
            ImGui::InputFloat(
                "Translation speed##Physics",
                &translation_speed,
                0.01f,
                0.1f,
                "%.2f");

            ImGui::Button("Translate##Physics");
            if (ImGui::IsItemActive())
            {
                translation += static_cast<sbs::scalar_type>(translation_speed);
                std::vector<Eigen::Vector3d> const dirichlet_boundary_conditions =
                    sbs::geometry::translate(
                        right_dirichlet_positions,
                        Eigen::Vector3d{translation, 0., 0.});

                simulation.apply_dirichlet_boundary_conditions(
                    right_bs,
                    right_fixed_particles,
                    dirichlet_boundary_conditions);
            }
            ImGui::Button("Untranslate##Physics");
            if (ImGui::IsItemActive())
            {
                translation -= static_cast<sbs::scalar_type>(translation_speed);
                std::vector<Eigen::Vector3d> const dirichlet_boundary_conditions =
                    sbs::geometry::translate(
                        right_dirichlet_positions,
                        Eigen::Vector3d{translation, 0., 0.});

                simulation.apply_dirichlet_boundary_conditions(
                    right_bs,
                    right_fixed_particles,
                    dirichlet_boundary_conditions);
            }
        }

        static std::vector<sbs::scalar_type> neighbour_distribution{};
        neighbour_distribution.clear();
        neighbour_distribution.reserve(meshless_model.point_count());

        if (ImPlot::BeginPlot("Neighbour Distribution##Plot"))
        {
            for (auto j = 0u; j < meshless_model.point_count(); ++j)
            {
                auto const N = mechanical_model.neighbours_of_meshless_node(j).size();
                neighbour_distribution.push_back(static_cast<sbs::scalar_type>(N));
            }
            ImPlot::PlotBars("N", neighbour_distribution.data(), neighbour_distribution.size());
            ImPlot::EndPlot();
        }

        ImGui::End();
    };

    renderer.on_pre_render = [&](sbs::physics::xpbd::simulation_t& s) {
        auto const& particles = s.particles()[beam_idx];
        if (should_render_points)
        {
            renderer.clear_points();
            for (auto j = 0; j < meshless_model.dof_count(); ++j)
            {
                auto const& p = particles[meshless_particle_index_offset + j];
                std::array<float, 3u> const vertex_position{
                    static_cast<float>(p.x().x()),
                    static_cast<float>(p.x().y()),
                    static_cast<float>(p.x().z())};
                std::array<float, 3u> const vertex_normal{0.f, 0.f, 0.f};
                std::array<float, 3u> const vertex_color{1.f, 0.f, 0.f};
                renderer.add_point(vertex_position, vertex_normal, vertex_color);
            }
        }
        if (should_render_active_fem_nodes)
        {
            renderer.clear_points();
            for (auto i = 0u; i < mechanical_model.dof_count(); ++i)
            {
                if (mechanical_model.has_basis_function(i))
                {
                    auto const& p = particles[fem_particle_index_offset + i];
                    std::array<float, 3u> const vertex_position{
                        static_cast<float>(p.x().x()),
                        static_cast<float>(p.x().y()),
                        static_cast<float>(p.x().z())};
                    std::array<float, 3u> const vertex_normal{0.f, 0.f, 0.f};
                    std::array<float, 3u> const vertex_color{1.f, 0.f, 0.f};
                    renderer.add_point(vertex_position, vertex_normal, vertex_color);
                }
            }
        }
        renderer.clear_lines();
        std::array<float, 3u> const p1{0.f, 0.f, 0.f};
        std::array<float, 3u> const n1{0.f, 0.f, 1.f};
        std::array<float, 3u> const c1{1.f, 0.f, 0.f};

        std::array<float, 3u> const p2{0.f, 1.f, 0.f};
        std::array<float, 3u> const n2{0.f, 0.f, 1.f};
        std::array<float, 3u> const c2{1.f, 0.f, 0.f};
        renderer.add_line(p1, n1, c1, p2, n2, c2);
    };

    renderer.launch();

    return 0;
}