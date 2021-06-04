#ifndef SBS_COMMON_NODE_H
#define SBS_COMMON_NODE_H

#include "mesh.h"

#include <string>

namespace sbs {
namespace common {

/**
 * @brief
 */
struct node_t
{
    std::string id;

    shared_vertex_mesh_t physical_model;
    shared_vertex_surface_mesh_t render_model;

    struct render_state_t
    {
        bool should_transfer_vertices = true;
        bool should_transfer_indices  = true;
        bool should_render_wireframe  = false;
    } render_state;

    unsigned int VBO;
    unsigned int VAO;
    unsigned int EBO;
};

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_NODE_H