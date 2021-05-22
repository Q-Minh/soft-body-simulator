#ifndef SBS_COMMON_MESH_H
#define SBS_COMMON_MESH_H

#include <vector>

namespace sbs {
namespace common {

struct shared_vertex_triangle_mesh_t
{
    struct position_t
    {
        double x, y, z;
    };

    struct triangle_t
    {
        std::uint32_t v1, v2, v3;
    };

    void rescale(
        position_t const& boxmin = position_t{-1., -1., -1.},
        position_t const& boxmax = position_t{+1., +1., +1.});

    std::vector<position_t> positions;
    std::vector<triangle_t> triangles;
};

struct shared_vertex_tetrahedral_mesh_t
{
    struct position_t
    {
        double x, y, z;
    };

    struct tetrahedron_t
    {
        std::uint32_t v1, v2, v3, v4;
    };

    void rescale(
        position_t const& boxmin = position_t{-1., -1., -1.},
        position_t const& boxmax = position_t{+1., +1., +1.});

    std::vector<position_t> positions;
    std::vector<tetrahedron_t> tetrahedra;

    shared_vertex_triangle_mesh_t extract_boundary_surface_mesh() const;
};

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_MESH_H