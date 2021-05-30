#include "io/load_scene.h"

namespace sbs {
namespace io {

common::scene_t load_scene(std::filesystem::path const& path)
{
    bool const should_continue = std::filesystem::exists(path) && path.has_filename() &&
                                 path.has_extension() && path.extension().string() == ".json";
    if (!should_continue)
    {
        return {};
    }

    nlohmann::json scene_specification{};

    std::ifstream ifs{path.string()};
    ifs >> scene_specification;

    common::scene_t scene{};

    auto const& lights_specification  = scene_specification["lights"];
    auto const& objects_specification = scene_specification["objects"];

    auto const& directional_light_specification = lights_specification["directional"];
    auto const& point_light_specification       = lights_specification["point"];

    /**
     * initialize directional light
     */
    {
        common::scene_t::directional_light_t light{};

        auto const& direction = directional_light_specification["direction"];
        light.dx              = direction["x"].get<float>();
        light.dy              = direction["y"].get<float>();
        light.dz              = direction["z"].get<float>();

        auto const& ambient  = directional_light_specification["ambient"];
        auto const& diffuse  = directional_light_specification["diffuse"];
        auto const& specular = directional_light_specification["specular"];

        light.ambient.r = ambient["r"].get<float>();
        light.ambient.g = ambient["g"].get<float>();
        light.ambient.b = ambient["b"].get<float>();

        light.diffuse.r = diffuse["r"].get<float>();
        light.diffuse.g = diffuse["g"].get<float>();
        light.diffuse.b = diffuse["b"].get<float>();

        light.specular.r = specular["r"].get<float>();
        light.specular.g = specular["g"].get<float>();
        light.specular.b = specular["b"].get<float>();

        light.specular.exp = specular["exp"].get<float>();

        scene.directional_light = light;
    }
    /**
     * initialize point light
     */
    {
        common::scene_t::point_light_t light{};

        auto const& position = point_light_specification["position"];
        light.x              = position["x"].get<float>();
        light.y              = position["y"].get<float>();
        light.z              = position["z"].get<float>();

        auto const& ambient  = point_light_specification["ambient"];
        auto const& diffuse  = point_light_specification["diffuse"];
        auto const& specular = point_light_specification["specular"];

        light.ambient.r = ambient["r"].get<float>();
        light.ambient.g = ambient["g"].get<float>();
        light.ambient.b = ambient["b"].get<float>();

        light.diffuse.r = diffuse["r"].get<float>();
        light.diffuse.g = diffuse["g"].get<float>();
        light.diffuse.b = diffuse["b"].get<float>();

        light.specular.r = specular["r"].get<float>();
        light.specular.g = specular["g"].get<float>();
        light.specular.b = specular["b"].get<float>();

        light.specular.exp = specular["exp"].get<float>();

        auto const& attenuation     = point_light_specification["attenuation"];
        light.attenuation.constant  = attenuation["constant"].get<float>();
        light.attenuation.linear    = attenuation["linear"].get<float>();
        light.attenuation.quadratic = attenuation["quadratic"].get<float>();

        scene.point_light = light;
    }

    for (auto const& object_spec : objects_specification)
    {
        auto const& geometry_spec           = object_spec["geometry"];
        std::string const geometry_type     = geometry_spec["type"].get<std::string>();
        std::string const asset_path_string = geometry_spec["path"].get<std::string>();

        std::filesystem::path const asset_path = path.parent_path() / asset_path_string;

        bool const can_process_file_format =
            asset_path.has_extension() && (asset_path.extension() == ".ply");

        bool const is_valid_file = std::filesystem::exists(asset_path) && asset_path.has_filename();

        bool const should_process_asset = is_valid_file && can_process_file_format;

        if (!should_process_asset)
            continue;

        auto const& physics_spec              = object_spec["physics"];
        std::string const body_type_str       = physics_spec["type"].get<std::string>();
        common::node_t::body_type_t body_type = body_type_str == "soft" ?
                                                    common::node_t::body_type_t::soft :
                                                    common::node_t::body_type_t::rigid;
        double const mass_per_vertex = physics_spec["mass"].get<double>();
        auto const& velocity_spec    = physics_spec["velocity"];
        double const vx              = velocity_spec["x"].get<double>();
        double const vy              = velocity_spec["y"].get<double>();
        double const vz              = velocity_spec["z"].get<double>();
        bool const is_fixed =
            physics_spec.contains("fixed") ? physics_spec["fixed"].get<bool>() : false;

        auto mesh_node          = std::make_shared<common::node_t>();
        mesh_node->id           = object_spec["id"].get<std::string>();
        mesh_node->body_type    = body_type;
        mesh_node->is_fixed     = is_fixed;
        mesh_node->render_state = common::node_t::render_state_t::dirty;

        std::optional<common::geometry_t> geometry = io::read_ply(asset_path);
        if (!geometry.has_value())
            continue;

        if (geometry->colors.empty())
        {
            geometry->colors.reserve(geometry->positions.size());
            auto const& color    = object_spec["color"];
            std::uint8_t const r = color["r"].get<std::uint8_t>();
            std::uint8_t const g = color["g"].get<std::uint8_t>();
            std::uint8_t const b = color["b"].get<std::uint8_t>();

            for (std::size_t i = 0u; i < geometry->positions.size(); i += 3u)
            {
                geometry->colors.push_back(r);
                geometry->colors.push_back(g);
                geometry->colors.push_back(b);
            }
        }

        common::shared_vertex_mesh_t mesh{geometry.value()};

        mesh.masses().setConstant(mass_per_vertex);

        mesh.velocities().row(0u).setConstant(vx);
        mesh.velocities().row(1u).setConstant(vy);
        mesh.velocities().row(2u).setConstant(vz);

        if (object_spec.contains("box"))
        {
            auto const& box_spec     = object_spec["box"];
            auto const& box_min_spec = box_spec["min"];
            auto const& box_max_spec = box_spec["max"];

            double const xmin = box_min_spec["x"].get<double>();
            double const ymin = box_min_spec["y"].get<double>();
            double const zmin = box_min_spec["z"].get<double>();
            double const xmax = box_max_spec["x"].get<double>();
            double const ymax = box_max_spec["y"].get<double>();
            double const zmax = box_max_spec["z"].get<double>();
            mesh.rescale({xmin, ymin, zmin}, {xmax, ymax, zmax});
        }
        if (object_spec.contains("translation"))
        {
            auto const& translation_spec = object_spec["translation"];
            double const tx              = translation_spec["x"].get<double>();
            double const ty              = translation_spec["y"].get<double>();
            double const tz              = translation_spec["z"].get<double>();

            mesh.positions().colwise() += Eigen::Vector3d{tx, ty, tz};
        }

        mesh_node->mesh = mesh;
        scene.objects.push_back(mesh_node);
    }

    return scene;
}

} // namespace io
} // namespace sbs