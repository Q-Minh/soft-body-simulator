#include "physics/xpbd/solver.h"
#include "rendering/renderer.h"

#include <chrono>
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