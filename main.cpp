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
    renderer.on_scene_loaded = [&](sbs::common::scene_t& scene) {
        std::vector<sbs::physics::xpbd::simulation_parameters_t> per_body_simulation_parameters{};
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

        solver.setup(scene.objects, per_body_simulation_parameters);
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
        for (auto const& body : scene.objects)
        {
            if (body->is_fixed)
                continue;

            /**
             * Reset external forces
             */
            body->mesh.forces().setZero();
            body->mesh.forces().colwise() += Eigen::Vector3d{0., -9.81, 0.};
        }

        solver.step(timestep, iterations, substeps);

        for (auto const& body : scene.objects)
        {
            if (body->is_fixed)
                continue;

            if (body->body_type == sbs::common::node_t::body_type_t::soft)
            {
                body->mesh.extract_boundary_surface_mesh();
                body->mesh.extract_boundary_normals();
            }
        }

        auto const end = std::chrono::steady_clock::now();
        auto const duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
    };

    renderer.on_new_imgui_frame = [](sbs::common::scene_t& scene) {
        ImGui::Begin("Soft Body Simulator");

        if (ImGui::CollapsingHeader("Physics", ImGuiTreeNodeFlags_DefaultOpen))
        {
            ImGui::TreePush();
            if (ImGui::CollapsingHeader("XPBD"))
            {
            }
            ImGui::TreePop();
        }
        if (ImGui::CollapsingHeader("Cutting", ImGuiTreeNodeFlags_DefaultOpen))
        {
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
                }
            }
            ImGui::TreePop();
        }

        ImGui::End();
    };

    bool const initialization_success = renderer.initialize(scene_specification_path);
    bool const shader_loading_success =
        renderer.use_shaders(vertex_shader_path, fragment_shader_path);

    if (initialization_success && shader_loading_success)
    {
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