#include <algorithm>
#include <array>
#include <imgui/imgui.h>
#include <implot/implot.h>
#include <iostream>
#include <sbs/geometry/get_simple_bar_model.h>
#include <sbs/geometry/get_simple_plane_model.h>
#include <sbs/math/kernels.h>
#include <sbs/math/mapping.h>
#include <sbs/math/quadrature.h>
#include <sbs/physics/body/environment_body.h>
#include <sbs/physics/body/fem_mixed_body.h>
#include <sbs/physics/collision/brute_force_cd_system.h>
#include <sbs/physics/mechanics/fem_sph_model.h>
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
using fem_mixed_model_type = sbs::physics::mechanics::fem_sph_model_t<kernel_function_type>;
using fem_model_type       = typename fem_mixed_model_type::fem_model_type;
using meshless_model_type  = typename fem_mixed_model_type::meshless_model_type;
using visual_model_type = sbs::physics::visual::fem_mixed_embedded_surface_t<fem_mixed_model_type>;
using body_type         = sbs::physics::body::fem_mixed_body_t<fem_mixed_model_type>;

sbs::scalar_type initial_volume(fem_mixed_model_type& body)
{
    sbs::scalar_type V{0.};
    auto const& topology = body.domain().topology();
    for (auto t = 0u; t < topology.tetrahedron_count(); ++t)
    {
        auto const v1 = topology.tetrahedron(t).v1();
        auto const v2 = topology.tetrahedron(t).v2();
        auto const v3 = topology.tetrahedron(t).v3();
        auto const v4 = topology.tetrahedron(t).v4();

        Eigen::Vector3d const& X1 = body.point(v1);
        Eigen::Vector3d const& X2 = body.point(v2);
        Eigen::Vector3d const& X3 = body.point(v3);
        Eigen::Vector3d const& X4 = body.point(v4);

        sbs::scalar_type const Vt = sbs::geometry::tetrahedron_volume(X1, X2, X3, X4);
        V += Vt;
    }
    return V;
}

sbs::scalar_type dof_volume(fem_mixed_model_type& body)
{
    sbs::scalar_type V{0.};
    auto const& topology = body.domain().topology();
    for (auto t = 0u; t < topology.tetrahedron_count(); ++t)
    {
        auto const v1 = topology.tetrahedron(t).v1();
        auto const v2 = topology.tetrahedron(t).v2();
        auto const v3 = topology.tetrahedron(t).v3();
        auto const v4 = topology.tetrahedron(t).v4();

        Eigen::Vector3d const& x1 = body.dof(v1);
        Eigen::Vector3d const& x2 = body.dof(v2);
        Eigen::Vector3d const& x3 = body.dof(v3);
        Eigen::Vector3d const& x4 = body.dof(v4);

        sbs::scalar_type const Vt = sbs::geometry::tetrahedron_volume(x1, x2, x3, x4);
        V += Vt;
    }
    return V;
}

sbs::scalar_type interpolated_volume(fem_mixed_model_type& body)
{
    sbs::scalar_type V{0.};
    auto const& topology = body.domain().topology();
    for (auto t = 0u; t < topology.tetrahedron_count(); ++t)
    {
        auto const v1 = topology.tetrahedron(t).v1();
        auto const v2 = topology.tetrahedron(t).v2();
        auto const v3 = topology.tetrahedron(t).v3();
        auto const v4 = topology.tetrahedron(t).v4();

        Eigen::Vector3d const& X1 = body.point(v1);
        Eigen::Vector3d const& X2 = body.point(v2);
        Eigen::Vector3d const& X3 = body.point(v3);
        Eigen::Vector3d const& X4 = body.point(v4);

        Eigen::Vector3d const& x1 = body.mixed_interpolation_field_at(X1).eval();
        Eigen::Vector3d const& x2 = body.mixed_interpolation_field_at(X2).eval();
        Eigen::Vector3d const& x3 = body.mixed_interpolation_field_at(X3).eval();
        Eigen::Vector3d const& x4 = body.mixed_interpolation_field_at(X4).eval();

        sbs::scalar_type const Vt = sbs::geometry::tetrahedron_volume(x1, x2, x3, x4);
        V += Vt;
    }
    return V;
}

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
    simulation.simulation_parameters().compliance                  = 1e-12;
    simulation.simulation_parameters().damping                     = 1e-6;
    simulation.simulation_parameters().collision_compliance        = 1e-5;
    simulation.simulation_parameters().collision_damping           = 1e-2;
    simulation.simulation_parameters().poisson_ratio               = 0.45;
    simulation.simulation_parameters().young_modulus               = 1e6;
    simulation.simulation_parameters().positional_penalty_strength = 4.;

    // Load geometry
    sbs::common::geometry_t beam_geometry = sbs::geometry::get_simple_bar_model(24u, 24u, 24u);
    beam_geometry.set_color(255, 255, 0);
    Eigen::Affine3d beam_transform{Eigen::Translation3d(-1., 4., 2.)};
    // beam_transform.rotate(
    //     Eigen::AngleAxisd(3.14159 / 20., Eigen::Vector3d{0., 1., 0.2}.normalized()));
    beam_transform.scale(Eigen::Vector3d{.4, 0.2, .4});
    beam_geometry                  = sbs::common::transform(beam_geometry, beam_transform);
    sbs::scalar_type const support = 2.;
    std::array<unsigned int, 3u> const resolution{24u, 24u, 24u};

    // Initialize soft body
    auto const beam_idx = simulation.add_body();
    simulation.bodies()[beam_idx] =
        std::make_unique<body_type>(simulation, beam_idx, beam_geometry, resolution, support);
    body_type& beam = *dynamic_cast<body_type*>(simulation.bodies()[beam_idx].get());

    fem_mixed_model_type& mechanical_model = beam.get_mechanical_model();
    sbs::scalar_type const mass_density    = 3.;

    // Create FEM particles
    auto const fem_particle_index_offset = simulation.particles()[beam_idx].size();
    for (auto i = 0u; i < mechanical_model.dof_count(); ++i)
    {
        Eigen::Vector3d const& xi = mechanical_model.dof(i);
        sbs::physics::xpbd::particle_t p{xi};
        // p.mass() = 0.;
        p.mass() = 1.;
        simulation.add_particle(p, beam_idx);
    }

    // Create SPH particles
    auto const meshless_particle_index_offset = simulation.particles()[beam_idx].size();
    auto& meshless_model                      = mechanical_model.meshless_model();
    for (auto j = 0u; j < meshless_model.dof_count(); ++j)
    {
        Eigen::Vector3d const& xj = meshless_model.dof(j);
        sbs::physics::xpbd::particle_t p{xj};
        // p.mass() = 0.;
        p.mass() = 1.;
        simulation.add_particle(p, beam_idx);
    }
    beam.set_fem_particle_offset(fem_particle_index_offset);
    beam.set_meshless_particle_offset(meshless_particle_index_offset);

    beam.get_visual_model().update();

    // Create mass distribution
    // auto const& domain   = mechanical_model.domain();
    // auto const& topology = domain.topology();
    // for (auto t = 0u; t < topology.tetrahedron_count(); ++t)
    //{
    //    std::vector<sbs::index_type> const& meshless_particles_in_tet =
    //        mechanical_model.particles_in_tetrahedron(t);
    //    auto N                                  = meshless_particles_in_tet.size();
    //    sbs::topology::tetrahedron_t const& tet = topology.tetrahedron(t);
    //    N += tet.vertex_indices().size();

    //    sbs::scalar_type const det =
    //        static_cast<sbs::scalar_type>(domain.barycentric_map(t).determinant());
    //    sbs::scalar_type const Vtet = (1. / 6.) * det;
    //    sbs::scalar_type const M    = (Vtet * mass_density);
    //    sbs::scalar_type const dM   = M / static_cast<sbs::scalar_type>(N);

    //    std::vector<sbs::physics::xpbd::particle_t>& beam_particles =
    //        simulation.particles()[beam_idx];
    //    for (sbs::index_type const i : tet.vertex_indices())
    //    {
    //        sbs::physics::xpbd::particle_t& p = beam_particles[fem_particle_index_offset + i];
    //        p.mass() += dM;
    //    }
    //    for (sbs::index_type const j : meshless_particles_in_tet)
    //    {
    //        sbs::physics::xpbd::particle_t& p = beam_particles[meshless_particle_index_offset +
    //        j]; p.mass() += dM;
    //    }
    //}

    auto const alpha = simulation.simulation_parameters().compliance;
    auto const beta  = simulation.simulation_parameters().damping;
    auto const nu    = simulation.simulation_parameters().poisson_ratio;
    auto const E     = simulation.simulation_parameters().young_modulus;

    // Create constraints
    auto const& domain   = mechanical_model.domain();
    auto const& topology = domain.topology();
    bool is_fully_mixed{true};
    for (auto e = 0u; e < mechanical_model.element_count(); ++e)
    {
        bool const is_mixed_tet = !mechanical_model.particles_in_tetrahedron(e).empty();
        if (is_mixed_tet)
        {
            using constraint_type =
                sbs::physics::xpbd::stvk_fem_mixed_nodal_integration_constraint_t<
                    fem_mixed_model_type>;

            std::vector<sbs::index_type> const& meshless_particles_in_tet =
                mechanical_model.particles_in_tetrahedron(e);
            auto const N = meshless_particles_in_tet.size();

            sbs::scalar_type const det =
                static_cast<sbs::scalar_type>(domain.barycentric_map(e).determinant());
            sbs::scalar_type const Vtet = (1. / 6.) * det;
            sbs::scalar_type const Vj   = Vtet / static_cast<sbs::scalar_type>(N);

            for (sbs::index_type const j : meshless_particles_in_tet)
            {
                auto constraint = std::make_unique<constraint_type>(
                    alpha,
                    beta,
                    j,
                    beam_idx,
                    Vj,
                    e,
                    mechanical_model,
                    E,
                    nu,
                    fem_particle_index_offset,
                    meshless_particle_index_offset);

                simulation.add_constraint(std::move(constraint));
            }
        }
        else
        {
            is_fully_mixed   = false;
            auto const& cell = mechanical_model.cell(e);

            auto const v1 = cell.node(0u);
            auto const v2 = cell.node(1u);
            auto const v3 = cell.node(2u);
            auto const v4 = cell.node(3u);

            auto const& element = mechanical_model.element(e);
            auto const& mapping = element.mapping();
            sbs::math::tetrahedron_1point_constant_quadrature_rule_t<decltype(mapping)>
                quadrature_rule{mapping};

            auto const& Xg = quadrature_rule.points.front().cast<sbs::scalar_type>();
            auto const& wg = static_cast<sbs::scalar_type>(quadrature_rule.weights.front());

            using constraint_type =
                sbs::physics::xpbd::stvk_tetrahedral_quadrature_strain_constraint_t<
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
    }

    if (is_fully_mixed)
        std::cout << "Full hybrid model\n";
    else
        std::cout << "Partial hybrid model\n";

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
    timestep.iterations() = 10u;
    timestep.substeps()   = 3u;
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
             button == GLFW_MOUSE_BUTTON_LEFT && mods == GLFW_MOD_CONTROL && action ==
             GLFW_PRESS);
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
                "FEM dofs: " + std::to_string(mechanical_model.dof_count());
            std::string const sph_dofs_str =
                "SPH dofs: " + std::to_string(meshless_model.dof_count());
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

            // sbs::scalar_type const V     = interpolated_volume(mechanical_model);
            // std::string const volume_str = "Volume: " + std::to_string(V);
            // ImGui::Text(volume_str.c_str());

            // sbs::scalar_type const V0         = initial_volume(mechanical_model);
            // std::string const rest_volume_str = "Rest volume: " + std::to_string(V0);
            // ImGui::Text(rest_volume_str.c_str());

            // sbs::scalar_type const Vd        = dof_volume(mechanical_model);
            // std::string const dof_volume_str = "Dof volume: " + std::to_string(Vd);
            // ImGui::Text(dof_volume_str.c_str());
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
        if (should_render_active_fem_nodes)
        {
            renderer.clear_points();
            for (auto i = 0u; i < mechanical_model.dof_count(); ++i)
            {
                if (mechanical_model.has_basis_function(i))
                {
                    auto const& p = particles[fem_particle_index_offset + i];
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