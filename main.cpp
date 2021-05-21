#include "rendering/renderer.h"

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
    bool const initialization_success = renderer.initialize(scene_specification_path);
    bool const shader_loading_success =
        renderer.use_shaders(vertex_shader_path, fragment_shader_path);

    if (initialization_success && shader_loading_success)
        renderer.launch();

    return 0;
}