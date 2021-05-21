#ifndef SBS_PHYSICS_NODE_H
#define SBS_PHYSICS_NODE_H

#include "common/node.h"
#include "mesh.h"

namespace sbs {
namespace physics {

struct node_t : common::node_t
{
    enum class body_type_t { soft, rigid };
    body_type_t body_type;
};

struct triangle_mesh_node_t : node_t
{
    virtual void prepare_for_rendering() override;

    physics::shared_vertex_triangle_mesh_t mesh;
};

struct tetrahedral_mesh_node_t : node_t
{
    virtual void prepare_for_rendering() override;

    physics::shared_vertex_tetrahedral_mesh_t mesh;
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_NODE_H