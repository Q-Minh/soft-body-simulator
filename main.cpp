#include "physics/cutting/cut_tetrahedron.h"
#include "physics/xpbd/solver.h"
#include "rendering/renderer.h"

#include <chrono>
#include <imgui/imgui.h>
#include <iostream>

int main(int argc, char** argv)
{
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

    /**
     * physics update goes here
     */
    std::vector<sbs::physics::xpbd::simulation_parameters_t> per_body_simulation_parameters{};
    renderer.on_scene_loaded = [&](sbs::common::scene_t& scene) {
        per_body_simulation_parameters.clear();
        per_body_simulation_parameters.reserve(scene.objects.size());

        for (auto const& body : scene.objects)
        {
            sbs::physics::xpbd::simulation_parameters_t params{};
            params.alpha           = 1e-3;
            params.constraint_type = sbs::physics::xpbd::constraint_type_t::distance;

            per_body_simulation_parameters.push_back(params);

            if (body->body_type == sbs::common::node_t::body_type_t::soft)
            {
                body->mesh.extract_boundary_surface_mesh();
            }
            body->mesh.extract_boundary_normals();
        }

        solver.setup(&scene.objects, per_body_simulation_parameters);
    };

    renderer.on_new_physics_timestep = [&](double render_frame_dt, sbs::common::scene_t& scene) {
        static double tb                   = 0.;
        double constexpr timestep          = 1. / 60.;
        std::uint32_t constexpr iterations = 60u;
        std::uint32_t constexpr substeps   = 60u;

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
        // for (auto const& body : scene.objects)
        //{
        //    if (body->is_fixed)
        //        continue;

        //    /**
        //     * Reset external forces
        //     */
        //    body->mesh.forces().setZero();
        //    body->mesh.forces().colwise() += Eigen::Vector3d{0., -9.81, 0.};
        //}

        // solver.step(timestep, iterations, substeps);

        for (auto const& body : scene.objects)
        {
            if (body->is_fixed)
                continue;

            if (body->body_type == sbs::common::node_t::body_type_t::soft)
            {
                body->mesh.extract_boundary_surface_mesh();
                body->mesh.extract_boundary_normals();
                body->render_state = sbs::common::node_t::render_state_t::dirty;
            }
        }

        auto const end = std::chrono::steady_clock::now();
        auto const duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    };

    renderer.on_new_imgui_frame = [&](sbs::common::scene_t& scene) {
        ImGui::Begin("Soft Body Simulator");

        static int active_body = 0;
        if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_None))
        {
            float const w = ImGui::GetColumnWidth();

            ImGui::BulletText("Select active object");
            for (std::size_t b = 0u; b < scene.objects.size(); ++b)
            {
                auto const& body   = scene.objects[b];
                auto const& params = per_body_simulation_parameters[b];

                ImGui::RadioButton(body->id.c_str(), &active_body, b);
            }

            if (ImGui::Button("Reload", ImVec2(w / 2.f, 0.f)))
            {
                renderer.load_scene(scene_specification_path);
            }
        }

        if (ImGui::CollapsingHeader("Physics", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::TreePush();
            if (ImGui::CollapsingHeader("XPBD"))
            {
            }
            ImGui::TreePop();
        }
        auto const& tet_body = scene.objects[active_body];
        bool const is_tet =
            (tet_body->mesh.elements().rows() == 4) && (tet_body->mesh.elements().cols() == 1);
        if (ImGui::CollapsingHeader("Cutting", ImGuiTreeNodeFlags_DefaultOpen) && is_tet)
        {
            auto const& P = tet_body->mesh.positions();
            auto const v1 = tet_body->mesh.elements()(0u, 0u);
            auto const v2 = tet_body->mesh.elements()(1u, 0u);
            auto const v3 = tet_body->mesh.elements()(2u, 0u);
            auto const v4 = tet_body->mesh.elements()(3u, 0u);

            auto const& p1 = P.col(v1);
            auto const& p2 = P.col(v2);
            auto const& p3 = P.col(v3);
            auto const& p4 = P.col(v4);

            std::array<Eigen::Vector3d, 6u> edge_intersections{
                0.5 * (p1 + p2),
                0.5 * (p2 + p3),
                0.5 * (p3 + p1),
                0.5 * (p1 + p4),
                0.5 * (p2 + p4),
                0.5 * (p3 + p4),
            };
            std::array<Eigen::Vector3d, 4u> face_intersections{
                0.33 * (p1 + p2 + p4),
                0.33 * (p2 + p3 + p4),
                0.33 * (p3 + p1 + p4),
                0.33 * (p1 + p3 + p2)};

            auto const update_scene_after_cut =
                [&renderer, &scene](
                    int cut_body,
                    std::vector<std::pair<
                        Eigen::Matrix3Xd,
                        Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic>>> const& tets) {
                    renderer.remove_object_from_scene(active_body);
                    for (auto const& [P, T] : tets)
                    {
                        sbs::common::shared_vertex_mesh_t mesh{P, T};
                        mesh.set_color({1.f, 1.f, 0.f});
                        auto const node    = std::make_shared<sbs::common::node_t>();
                        node->body_type    = sbs::common::node_t::body_type_t::soft;
                        node->id           = "cut tet";
                        node->is_fixed     = false;
                        node->render_state = sbs::common::node_t::render_state_t::dirty;
                        node->mesh         = std::move(mesh);
                        scene.objects.push_back(node);
                    }
                    renderer.load_scene(scene);
                };

            ImGui::TreePush();
            if (ImGui::CollapsingHeader("Common Case 1"))
            {
                float const w     = ImGui::GetColumnWidth();
                static int choice = 0;
                ImGui::RadioButton("1,2,5##Case1", &choice, 0);
                ImGui::RadioButton("1,3,4##Case1", &choice, 1);
                ImGui::RadioButton("2,3,6##Case1", &choice, 2);
                ImGui::RadioButton("4,5,6##Case1", &choice, 3);
                if (ImGui::Button("Cut##Case1", ImVec2(w / 2.f, 0.f)))
                {
                    std::vector<std::pair<
                        Eigen::Matrix3Xd,
                        Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic>>>
                        tets{};

                    if (choice == 0)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00010011},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 1)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00001101},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 2)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00100110},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 3)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00111000},
                            edge_intersections,
                            face_intersections);
                    }

                    update_scene_after_cut(active_body, tets);
                }
            }
            ImGui::TreePop();

            ImGui::TreePush();
            if (ImGui::CollapsingHeader("Common Case 2"))
            {
                float const w     = ImGui::GetColumnWidth();
                static int choice = 0;
                ImGui::RadioButton("1,2,4,6##Case2", &choice, 0);
                ImGui::RadioButton("1,3,5,6##Case2", &choice, 1);
                ImGui::RadioButton("2,3,4,5##Case2", &choice, 2);
                if (ImGui::Button("Cut##Case1", ImVec2(w / 2.f, 0.f)))
                {
                    std::vector<std::pair<
                        Eigen::Matrix3Xd,
                        Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic>>>
                        tets{};

                    if (choice == 0)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00101011},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 1)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00110101},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 2)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00011110},
                            edge_intersections,
                            face_intersections);
                    }

                    update_scene_after_cut(active_body, tets);
                }
            }
            ImGui::TreePop();

            ImGui::TreePush();
            if (ImGui::CollapsingHeader("Common Case 3"))
            {
                float const w     = ImGui::GetColumnWidth();
                static int choice = 0;
                ImGui::RadioButton("1##Case3", &choice, 0);
                ImGui::RadioButton("2##Case3", &choice, 1);
                ImGui::RadioButton("3##Case3", &choice, 2);
                ImGui::RadioButton("4##Case3", &choice, 3);
                ImGui::RadioButton("5##Case3", &choice, 4);
                ImGui::RadioButton("6##Case3", &choice, 5);
                if (ImGui::Button("Cut##Case1", ImVec2(w / 2.f, 0.f)))
                {
                    std::vector<std::pair<
                        Eigen::Matrix3Xd,
                        Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic>>>
                        tets{};

                    if (choice == 0)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00000001},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 1)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00000010},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 2)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00000100},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 3)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00001000},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 4)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00010000},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 5)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00100000},
                            edge_intersections,
                            face_intersections);
                    }

                    update_scene_after_cut(active_body, tets);
                }
            }
            ImGui::TreePop();

            ImGui::TreePush();
            if (ImGui::CollapsingHeader("Common Case 4"))
            {
                float const w     = ImGui::GetColumnWidth();
                static int choice = 0;
                ImGui::RadioButton("1,2##Case4", &choice, 0);
                ImGui::RadioButton("1,3##Case4", &choice, 1);
                ImGui::RadioButton("1,4##Case4", &choice, 2);
                ImGui::RadioButton("1,5##Case4", &choice, 3);
                ImGui::RadioButton("2,3##Case4", &choice, 4);
                ImGui::RadioButton("2,5##Case4", &choice, 5);
                ImGui::RadioButton("2,6##Case4", &choice, 6);
                ImGui::RadioButton("3,4##Case4", &choice, 7);
                ImGui::RadioButton("3,6##Case4", &choice, 8);
                ImGui::RadioButton("4,5##Case4", &choice, 9);
                ImGui::RadioButton("4,6##Case4", &choice, 10);
                ImGui::RadioButton("5,6##Case4", &choice, 11);
                if (ImGui::Button("Cut##Case1", ImVec2(w / 2.f, 0.f)))
                {
                    std::vector<std::pair<
                        Eigen::Matrix3Xd,
                        Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic>>>
                        tets{};

                    if (choice == 0)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00000011},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 1)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00000101},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 2)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00001001},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 3)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00010001},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 4)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00000110},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 5)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00010010},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 6)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00100010},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 7)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00001100},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 8)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00100100},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 9)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00011000},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 10)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00101000},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 11)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00110000},
                            edge_intersections,
                            face_intersections);
                    }

                    update_scene_after_cut(active_body, tets);
                }
            }
            ImGui::TreePop();

            ImGui::TreePush();
            if (ImGui::CollapsingHeader("Common Case 5"))
            {
                float const w     = ImGui::GetColumnWidth();
                static int choice = 0;
                ImGui::RadioButton("1,2,4##Case5", &choice, 0);
                ImGui::RadioButton("1,2,6##Case5", &choice, 1);
                ImGui::RadioButton("1,3,5##Case5", &choice, 2);
                ImGui::RadioButton("1,3,6##Case5", &choice, 3);
                ImGui::RadioButton("1,4,6##Case5", &choice, 4);
                ImGui::RadioButton("1,5,6##Case5", &choice, 5);
                ImGui::RadioButton("2,3,4##Case5", &choice, 6);
                ImGui::RadioButton("2,3,5##Case5", &choice, 7);
                ImGui::RadioButton("2,4,5##Case5", &choice, 8);
                ImGui::RadioButton("2,4,6##Case5", &choice, 9);
                ImGui::RadioButton("3,4,5##Case5", &choice, 10);
                ImGui::RadioButton("3,5,6##Case5", &choice, 11);
                if (ImGui::Button("Cut##Case1", ImVec2(w / 2.f, 0.f)))
                {
                    std::vector<std::pair<
                        Eigen::Matrix3Xd,
                        Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic>>>
                        tets{};

                    if (choice == 0)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00001011},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 1)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00100011},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 2)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00010101},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 3)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00100101},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 4)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00101001},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 5)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00110001},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 6)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00001110},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 7)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00010110},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 8)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00011010},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 9)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00101010},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 10)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00011100},
                            edge_intersections,
                            face_intersections);
                    }
                    else if (choice == 11)
                    {
                        tets = sbs::physics::cutting::cut_tetrahedron(
                            tet_body->mesh.positions(),
                            tet_body->mesh.elements(),
                            0u,
                            std::byte{0b00110100},
                            edge_intersections,
                            face_intersections);
                    }

                    update_scene_after_cut(active_body, tets);
                }
            }
            ImGui::TreePop();
        }

        ImGui::End();
    };

    bool const initialization_success = renderer.initialize();
    bool const shader_loading_success =
        renderer.use_shaders(vertex_shader_path, fragment_shader_path);

    if (initialization_success && shader_loading_success)
    {
        renderer.load_scene(scene_specification_path);
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