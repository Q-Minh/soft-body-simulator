#include "physics/node.h"

#include <glm/glm.hpp>
#include <numeric>
#include <unordered_map>

namespace sbs {
namespace physics {

void triangle_mesh_node_t::prepare_for_rendering()
{
    this->positions.clear();
    this->positions.reserve(this->mesh.positions.size() * 3u);
    for (auto const& position : this->mesh.positions)
    {
        this->positions.push_back(position.x);
        this->positions.push_back(position.y);
        this->positions.push_back(position.z);
    }

    bool const should_compute_normals = this->normals.empty();
    if (should_compute_normals)
    {
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
    }

    // assumes topology changes (cutting)
    this->indices.clear();
    this->indices.reserve(this->mesh.triangles.size() * 3u);
    for (auto const& triangle : this->mesh.triangles)
    {
        this->indices.push_back(triangle.v1);
        this->indices.push_back(triangle.v2);
        this->indices.push_back(triangle.v3);
    }
}

void tetrahedral_mesh_node_t::prepare_for_rendering() {}

} // namespace physics
} // namespace sbs