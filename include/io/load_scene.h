#ifndef SBS_IO_LOAD_SCENE_H
#define SBS_IO_LOAD_SCENE_H

#include "common/scene.h"
#include "io/ply.h"

#include <filesystem>
#include <fstream>
#include <nlohmann/json.hpp>

namespace sbs {
namespace io {

common::scene_t load_scene(std::filesystem::path const& path);

} // namespace io
} // namespace sbs

#endif // SBS_IO_LOAD_SCENE_H