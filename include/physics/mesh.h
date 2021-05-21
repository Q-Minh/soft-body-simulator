#ifndef SBS_PHYSICS_MESH_H
#define SBS_PHYSICS_MESH_H

#include "common/mesh.h"

namespace sbs {
namespace physics {

struct shared_vertex_triangle_mesh_t : public common::shared_vertex_triangle_mesh_t
{
    struct velocity_t
    {
        double vx, vy, vz;
    };

    std::vector<double> masses;
    std::vector<velocity_t> velocities;
};

struct shared_vertex_tetrahedral_mesh_t : public common::shared_vertex_tetrahedral_mesh_t
{
    struct velocity_t
    {
        double vx, vy, vz;
    };

    std::vector<double> masses;
    std::vector<velocity_t> velocities;
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MESH_H