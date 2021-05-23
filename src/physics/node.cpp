#include "physics/node.h"

#include <algorithm>
#include <array>
#include <glm/glm.hpp>
#include <numeric>
#include <unordered_map>

namespace sbs {
namespace physics {

static void from_triangle_mesh(
    common::shared_vertex_triangle_mesh_t const& mesh,
    std::vector<float>& positions,
    std::vector<float>& normals,
    std::vector<std::uint32_t>& indices)
{
    /**
     * We assume the mesh deforms and might change topology constantly, so we always rebuild
     * positions, normals and topology.
     */
    positions.clear();
    normals.clear();
    indices.clear();

    positions.reserve(mesh.positions.size() * 3u);
    for (auto const& position : mesh.positions)
    {
        positions.push_back(position.x);
        positions.push_back(position.y);
        positions.push_back(position.z);
    }

    std::unordered_map<std::uint32_t, std::vector<glm::vec3>> one_ring_neighbour_normals{};

    for (auto const& triangle : mesh.triangles)
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
        normals.push_back(vertex_normal.x);
        normals.push_back(vertex_normal.y);
        normals.push_back(vertex_normal.z);
    }

    indices.reserve(mesh.triangles.size() * 3u);
    for (auto const& triangle : mesh.triangles)
    {
        indices.push_back(triangle.v1);
        indices.push_back(triangle.v2);
        indices.push_back(triangle.v3);
    }
}

void triangle_mesh_node_t::prepare_for_rendering()
{
    from_triangle_mesh(this->mesh, this->positions, this->normals, this->indices);
}

void tetrahedral_mesh_node_t::prepare_for_rendering()
{
    auto const boundary_mesh = this->mesh.extract_boundary_surface_mesh();
    from_triangle_mesh(boundary_mesh, this->positions, this->normals, this->indices);
}

} // namespace physics
} // namespace sbs