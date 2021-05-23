#include "common/mesh.h"

#include <algorithm>
#include <array>
#include <map>
#include <numeric>

namespace sbs {
namespace common {

template <class PositionType>
static void rescale_internal(
    PositionType const& boxmin,
    PositionType const& boxmax,
    std::vector<PositionType>& positions)
{
    using position_t = PositionType;

    auto const min_reduce_op = [](position_t const& current_min, position_t const& p) {
        position_t new_min{current_min};
        if (p.x < current_min.x)
            new_min.x = p.x;
        if (p.y < current_min.y)
            new_min.y = p.y;
        if (p.z < current_min.z)
            new_min.z = p.z;

        return new_min;
    };

    auto const max_reduce_op = [](position_t const& current_max, position_t const& p) {
        position_t new_max{current_max};
        if (p.x > current_max.x)
            new_max.x = p.x;
        if (p.y > current_max.y)
            new_max.y = p.y;
        if (p.z > current_max.z)
            new_max.z = p.z;

        return new_max;
    };

    position_t const min_position =
        std::reduce(positions.begin(), positions.end(), positions.front(), min_reduce_op);

    position_t const max_position =
        std::reduce(positions.begin(), positions.end(), positions.front(), max_reduce_op);

    double const dx = max_position.x - min_position.x;
    double const dy = max_position.y - min_position.y;
    double const dz = max_position.z - min_position.z;

    double constexpr eps           = 1e-8;
    bool const dx_division_by_zero = std::abs(dx) < eps;
    bool const dy_division_by_zero = std::abs(dy) < eps;
    bool const dz_division_by_zero = std::abs(dz) < eps;

    auto const map_to_new_box = [=](position_t const& pos) {
        position_t rescaled_p{};
        rescaled_p.x = dx_division_by_zero ?
                           pos.x :
                           boxmin.x + (boxmax.x - boxmin.x) * (pos.x - min_position.x) / dx;
        rescaled_p.y = dy_division_by_zero ?
                           pos.y :
                           boxmin.y + (boxmax.y - boxmin.y) * (pos.y - min_position.y) / dy;
        rescaled_p.z = dz_division_by_zero ?
                           pos.z :
                           boxmin.z + (boxmax.z - boxmin.z) * (pos.z - min_position.z) / dz;
        return rescaled_p;
    };

    std::transform(
        positions.begin(),
        positions.end(),
        positions.begin(),
        [map_to_new_box](position_t const& pos) { return map_to_new_box(pos); });
}

void shared_vertex_triangle_mesh_t::rescale(position_t const& boxmin, position_t const& boxmax)
{
    rescale_internal(boxmin, boxmax, this->positions);
}

void shared_vertex_tetrahedral_mesh_t::rescale(position_t const& boxmin, position_t const& boxmax)
{
    rescale_internal(boxmin, boxmax, this->positions);
}

shared_vertex_triangle_mesh_t
shared_vertex_tetrahedral_mesh_t::extract_boundary_surface_mesh() const
{
    using triangle_type = std::array<std::uint32_t, 3u>;
    std::vector<triangle_type> triangles{};
    triangles.reserve(this->tetrahedra.size() * 4u);

    /**
     * Extract all triangles from tet mesh (all 4 faces of all tets)
     */
    for (auto const& tetrahedron : this->tetrahedra)
    {
        triangle_type f1{}, f2{}, f3{}, f4{};

        f1[0] = tetrahedron.v1;
        f1[1] = tetrahedron.v3;
        f1[2] = tetrahedron.v2;

        f2[0] = tetrahedron.v1;
        f2[1] = tetrahedron.v2;
        f2[2] = tetrahedron.v4;

        f3[0] = tetrahedron.v2;
        f3[1] = tetrahedron.v3;
        f3[2] = tetrahedron.v4;

        f4[0] = tetrahedron.v1;
        f4[1] = tetrahedron.v4;
        f4[2] = tetrahedron.v3;

        triangles.push_back(f1);
        triangles.push_back(f2);
        triangles.push_back(f3);
        triangles.push_back(f4);
    }

    /**
     * Reorder triangle indices in the same order for all triangles
     */
    std::vector<triangle_type> triangle_copies{triangles};
    for (auto& triangle : triangle_copies)
    {
        std::sort(triangle.begin(), triangle.end());
    }

    /**
     * Count number of times each triangle is present in the tet mesh.
     */
    std::map<triangle_type, std::uint32_t> triangle_occurrences{};
    for (auto const& triangle : triangle_copies)
    {
        if (triangle_occurrences.find(triangle) == triangle_occurrences.end())
        {
            triangle_occurrences[triangle] = 1u;
        }
        else
        {
            ++triangle_occurrences[triangle];
        }
    }

    shared_vertex_triangle_mesh_t boundary_mesh{};
    for (std::size_t i = 0u; i < triangles.size(); ++i)
    {
        /**
         * Boundary triangles should only be present once in the tet mesh, since
         * there is no adjacent face to it.
         */
        if (triangle_occurrences[triangle_copies[i]] != 1u)
            continue;

        boundary_mesh.triangles.push_back(shared_vertex_triangle_mesh_t::triangle_t{
            triangles[i][0],
            triangles[i][1],
            triangles[i][2]});
    }

    for (auto const& position : this->positions)
    {
        boundary_mesh.positions.push_back(
            shared_vertex_triangle_mesh_t::position_t{position.x, position.y, position.z});
    }

    return boundary_mesh;
}

} // namespace common
} // namespace sbs