#include <imgui/imgui.h>
#include <sbs/geometry/curves.h>
#include <sbs/geometry/line.h>
#include <sbs/physics/xpbd/simulation.h>
#include <sbs/rendering/renderer.h>

int main(int argc, char** argv)
{
    sbs::physics::xpbd::simulation_t simulation{};

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

    sbs::geometry::piecewise_bezier_curve_t<3, 2> parametric_curve{};

    sbs::geometry::bezier_curve_t<3, 2> bezier_curve_0{};
    bezier_curve_0.set_control_points(
        {Eigen::Vector3d{0., 0.5, 0.}, Eigen::Vector3d{1., 1., 0.}, Eigen::Vector3d{2., 0., 0.}});
    sbs::geometry::bezier_curve_t<3, 2> bezier_curve_1{};
    bezier_curve_1.set_control_points(
        {Eigen::Vector3d{2., 0., 0.}, Eigen::Vector3d{3., -1., 0.}, Eigen::Vector3d{4., 0., 0.}});

    parametric_curve.add_segment(bezier_curve_0);
    parametric_curve.add_segment(bezier_curve_1);

    sbs::geometry::line_segment_t<3> starting_line_segment{
        Eigen::Vector3d{0., 0.5, 0.},
        Eigen::Vector3d{0., 0., 1.5}};

    auto swept_surface =
        std::make_unique<sbs::geometry::swept_line_segment_t>(starting_line_segment);
    swept_surface->set_dt(0.1);
    swept_surface->set_color(Eigen::Vector3f{0.f, 1.f, 1.f});
    swept_surface->set_parametric_curve(
        [parametric_curve](sbs::scalar_type t) -> std::optional<Eigen::Vector3d> {
            return parametric_curve.eval(t);
        });

    renderer.add_rendered_object(std::move(swept_surface));
    auto swept_surface_ptr =
        dynamic_cast<sbs::geometry::swept_line_segment_t*>(renderer.rendered_object(0u).get());

    renderer.on_new_imgui_frame = [&](sbs::physics::xpbd::simulation_t& s) {
        ImGui::Begin("Swept line segment surface");

        static bool should_render_wireframe = false;
        ImGui::Checkbox("Wireframe", &should_render_wireframe);
        if (should_render_wireframe)
        {
            swept_surface_ptr->mark_should_render_wireframe();
        }
        else
        {
            swept_surface_ptr->mark_should_render_triangles();
        }

        static float dt = 0.1f;
        ImGui::InputFloat("dt", &dt, 0.01f, 0.1f, "%.2f");
        swept_surface_ptr->set_dt(static_cast<sbs::scalar_type>(dt));

        static float timestep = 0.5f;
        ImGui::InputFloat("timestep", &timestep, 0.01f, 0.1f, "%.2f");

        if (ImGui::Button("Sweep"))
        {
            swept_surface_ptr->sweep(timestep);
            swept_surface_ptr->mark_vertices_dirty();
            swept_surface_ptr->mark_indices_dirty();
        }

        ImGui::End();
    };

    renderer.launch();

    return 0;
}