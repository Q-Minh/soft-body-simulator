#ifndef SBS_IO_LOAD_SCENE_H
#define SBS_IO_LOAD_SCENE_H

#include "common/scene.h"
#include "common/geometry.h"

#include <filesystem>
#include <functional>
#include <string>

namespace sbs {
namespace io {
namespace scene {

struct scene_body_info
{
    std::string id;
    common::geometry_t geometry;
};

struct physics_body_info : scene_body_info
{
    double mass_density;
    struct velocity_t
    {
        double vx, vy, vz;
    } velocity;
};

} // namespace scene

common::scene_t load_scene(
    std::filesystem::path const& path,
    std::function<std::shared_ptr<common::renderable_node_t>(scene::scene_body_info const&)>
        environment_body_factory,
    std::function<std::shared_ptr<common::renderable_node_t>(scene::physics_body_info const&)>
        physics_body_factory);

} // namespace io
} // namespace sbs

#endif // SBS_IO_LOAD_SCENE_H