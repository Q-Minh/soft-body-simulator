#include "io/load_scene.h"

#include "io/ply.h"

#include <fstream>
#include <nlohmann/json.hpp>

namespace sbs {
namespace io {

static void rescale(
    common::geometry_t& geometry,
    double bxmin,
    double bymin,
    double bzmin,
    double bxmax,
    double bymax,
    double bzmax)
{
    double xmin, ymin, zmin;
    xmin = ymin = zmin = std::numeric_limits<float>::max();
    double xmax, ymax, zmax;
    xmax = ymax = zmax = std::numeric_limits<float>::lowest();

    for (std::size_t i = 0u; i < geometry.positions.size(); i += 3u)
    {
        float const x = geometry.positions[i];
        float const y = geometry.positions[i + 1u];
        float const z = geometry.positions[i + 2u];

        if (x < xmin)
            xmin = x;
        if (y < ymin)
            ymin = y;
        if (z < zmin)
            zmin = z;
        if (x > xmax)
            xmax = x;
        if (y > ymax)
            ymax = y;
        if (z > zmax)
            zmax = z;
    }

    double const dx = xmax - xmin;
    double const dy = ymax - ymin;
    double const dz = zmax - zmin;

    double constexpr eps           = 1e-8;
    bool const dx_division_by_zero = std::abs(dx) < eps;
    bool const dy_division_by_zero = std::abs(dy) < eps;
    bool const dz_division_by_zero = std::abs(dz) < eps;

    for (std::size_t i = 0u; i < geometry.positions.size(); i += 3u)
    {
        float const x0 = geometry.positions[i];
        float const y0 = geometry.positions[i + 1u];
        float const z0 = geometry.positions[i + 2u];

        auto const x = dx_division_by_zero ? x0 : bxmin + (bxmax - bxmin) * (x0 - xmin) / dx;
        auto const y = dy_division_by_zero ? y0 : bymin + (bymax - bymin) * (y0 - ymin) / dy;
        auto const z = dz_division_by_zero ? z0 : bzmin + (bzmax - bzmin) * (z0 - zmin) / dz;

        geometry.positions[i]      = static_cast<float>(x);
        geometry.positions[i + 1u] = static_cast<float>(y);
        geometry.positions[i + 2u] = static_cast<float>(z);
    }
}

static void translate(common::geometry_t& geometry, double tx, double ty, double tz)
{
    for (std::size_t i = 0u; i < geometry.positions.size(); i += 3u)
    {
        float const x              = geometry.positions[i];
        float const y              = geometry.positions[i + 1u];
        float const z              = geometry.positions[i + 2u];
        geometry.positions[i]      = x + static_cast<float>(tx);
        geometry.positions[i + 1u] = y + static_cast<float>(ty);
        geometry.positions[i + 2u] = z + static_cast<float>(tz);
    }
}

common::scene_t load_scene(
    std::filesystem::path const& path,
    std::function<std::shared_ptr<common::renderable_node_t>(scene::scene_body_info const&)>
        environment_body_factory,
    std::function<std::shared_ptr<common::renderable_node_t>(scene::physics_body_info const&)>
        physics_body_factory)
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

    auto const& lights_specification            = scene_specification["lights"];
    auto const& directional_light_specification = lights_specification["directional"];
    auto const& point_light_specification       = lights_specification["point"];

    auto const& environment_specification = scene_specification["environment"];
    auto const& objects_specification     = scene_specification["objects"];

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

    for (auto const& environment_object_spec : environment_specification)
    {
        auto const& geometry_spec           = environment_object_spec["geometry"];
        std::string const geometry_type     = geometry_spec["type"].get<std::string>();
        std::string const asset_path_string = geometry_spec["path"].get<std::string>();

        std::filesystem::path const asset_path = path.parent_path() / asset_path_string;

        bool const can_process_file_format =
            asset_path.has_extension() && (asset_path.extension() == ".ply");

        bool const is_valid_file = std::filesystem::exists(asset_path) && asset_path.has_filename();

        bool const should_process_asset = is_valid_file && can_process_file_format;

        if (!should_process_asset)
            continue;

        std::optional<common::geometry_t> geometry = io::read_ply(asset_path);
        if (!geometry.has_value())
            continue;

        if (!geometry->has_colors())
        {
            geometry->colors.reserve(geometry->positions.size());
            auto const& color    = environment_object_spec["color"];
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

        if (environment_object_spec.contains("box"))
        {
            auto const& box_spec     = environment_object_spec["box"];
            auto const& box_min_spec = box_spec["min"];
            auto const& box_max_spec = box_spec["max"];

            double const bxmin = box_min_spec["x"].get<double>();
            double const bymin = box_min_spec["y"].get<double>();
            double const bzmin = box_min_spec["z"].get<double>();
            double const bxmax = box_max_spec["x"].get<double>();
            double const bymax = box_max_spec["y"].get<double>();
            double const bzmax = box_max_spec["z"].get<double>();

            rescale(*geometry, bxmin, bymin, bzmin, bxmax, bymax, bzmax);
        }
        if (environment_object_spec.contains("translation"))
        {
            auto const& translation_spec = environment_object_spec["translation"];
            double const tx              = translation_spec["x"].get<double>();
            double const ty              = translation_spec["y"].get<double>();
            double const tz              = translation_spec["z"].get<double>();

            translate(*geometry, tx, ty, tz);
        }

        std::string const id = environment_object_spec["id"].get<std::string>();

        scene::scene_body_info sbi;
        sbi.id       = id;
        sbi.geometry = geometry.value();
        auto node    = environment_body_factory(sbi);
        node->set_id(sbi.id);
        node->set_as_environment_body();
        scene.nodes.push_back(node);
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

        auto const& physics_spec        = object_spec["physics"];
        std::string const body_type_str = physics_spec["type"].get<std::string>();
        double const mass_density       = physics_spec["mass"].get<double>();
        auto const& velocity_spec       = physics_spec["velocity"];
        double const vx                 = velocity_spec["x"].get<double>();
        double const vy                 = velocity_spec["y"].get<double>();
        double const vz                 = velocity_spec["z"].get<double>();

        std::optional<common::geometry_t> geometry = io::read_ply(asset_path);
        if (!geometry.has_value())
            continue;

        std::string const id = object_spec["id"].get<std::string>();

        if (!geometry->has_colors())
        {
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

            rescale(*geometry, xmin, ymin, zmin, xmax, ymax, zmax);
        }
        if (object_spec.contains("translation"))
        {
            auto const& translation_spec = object_spec["translation"];
            double const tx              = translation_spec["x"].get<double>();
            double const ty              = translation_spec["y"].get<double>();
            double const tz              = translation_spec["z"].get<double>();

            translate(*geometry, tx, ty, tz);
        }

        scene::physics_body_info pbi;
        pbi.id           = id;
        pbi.geometry     = geometry.value();
        pbi.velocity.vx  = vx;
        pbi.velocity.vy  = vy;
        pbi.velocity.vz  = vz;
        pbi.mass_density = mass_density;

        auto node = physics_body_factory(pbi);
        node->set_id(pbi.id);
        node->set_as_physically_simulated_body();
        scene.nodes.push_back(node);
    }

    return scene;
}

} // namespace io
} // namespace sbs