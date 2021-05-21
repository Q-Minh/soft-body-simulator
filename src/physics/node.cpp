#include "physics/node.h"

#include <algorithm>
#include <array>
#include <glm/glm.hpp>
#include <numeric>
#include <unordered_map>

namespace sbs {
namespace physics {

void triangle_mesh_node_t::prepare_for_rendering()
{
    /**
     * We assume the mesh deforms and might change topology constantly, so we always rebuild
     * positions, normals and topology.
     */
    this->positions.clear();
    this->normals.clear();
    this->indices.clear();

    this->positions.reserve(this->mesh.positions.size() * 3u);
    for (auto const& position : this->mesh.positions)
    {
        this->positions.push_back(position.x);
        this->positions.push_back(position.y);
        this->positions.push_back(position.z);
    }

    std::unordered_map<std::uint32_t, std::vector<glm::vec3>> one_ring_neighbour_normals{};

    for (auto const& triangle : this->mesh.triangles)
    {
        auto const& v1 = triangle.v1;
        auto const& v2 = triangle.v2;
        auto const& v3 = triangle.v3;

        glm::vec3 const v21{
            mesh.positions[v2].x - mesh.positions[v1].x,
            mesh.positions[v2].y - mesh.positions[v1].y,
            mesh.positions[v2].z - mesh.positions[v1].z};

        glm::vec3 const v31{
            mesh.positions[v3].x - mesh.positions[v1].x,
            mesh.positions[v3].y - mesh.positions[v1].y,
            mesh.positions[v3].z - mesh.positions[v1].z};

        glm::vec3 const triangle_normal = glm::normalize(glm::cross(v21, v31));

        one_ring_neighbour_normals[v1].push_back(triangle_normal);
        one_ring_neighbour_normals[v2].push_back(triangle_normal);
        one_ring_neighbour_normals[v3].push_back(triangle_normal);
    }

    for (std::size_t vi = 0u; vi < mesh.positions.size(); ++vi)
    {
        auto const& position                            = mesh.positions[vi];
        std::vector<glm::vec3> const& neighbour_normals = one_ring_neighbour_normals[vi];
        auto const sum = std::reduce(neighbour_normals.begin(), neighbour_normals.end());
        glm::vec3 const vertex_normal = glm::normalize(sum);
        this->normals.push_back(vertex_normal.x);
        this->normals.push_back(vertex_normal.y);
        this->normals.push_back(vertex_normal.z);
    }

    this->indices.reserve(this->mesh.triangles.size() * 3u);
    for (auto const& triangle : this->mesh.triangles)
    {
        this->indices.push_back(triangle.v1);
        this->indices.push_back(triangle.v2);
        this->indices.push_back(triangle.v3);
    }
}

void tetrahedral_mesh_node_t::prepare_for_rendering()
{
    this->positions.clear();
    this->normals.clear();
    this->indices.clear();

    using triangle_type = std::array<std::uint32_t, 3u>;
    std::vector<triangle_type> triangles{};
    triangles.reserve(this->mesh.tetrahedra.size() * 4u);

    for (auto const& tetrahedron : this->mesh.tetrahedra)
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

    auto const rotate_triangle_indices = [](triangle_type const& f) {
        return triangle_type{f[1], f[2], f[0]};
    };

    auto const is_same_triangle =
        [rotate_triangle_indices](triangle_type const& f1, triangle_type const& f2) {
            triangle_type shifted_triangle{f2};
            if (f1 == shifted_triangle)
                return true;

            shifted_triangle = rotate_triangle_indices(shifted_triangle);
            if (f1 == shifted_triangle)
                return true;

            shifted_triangle = rotate_triangle_indices(shifted_triangle);
            if (f1 == shifted_triangle)
                return true;

            return false;
        };


}

} // namespace physics
} // namespace sbs