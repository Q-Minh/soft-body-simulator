#ifndef SBS_COMMON_NODE_H
#define SBS_COMMON_NODE_H

#include <vector>
#include "mesh.h"

namespace sbs {
namespace common {

/**
 * @brief
 */
struct node_t
{
    enum class body_type_t { soft, rigid };
    body_type_t body_type;

    shared_vertex_mesh_t mesh;

    enum class render_state_t { dirty, clean };
    render_state_t render_state;

    unsigned int VBO;
    unsigned int VAO;
    unsigned int EBO;
};

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_NODE_H