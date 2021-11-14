#include <Eigen/Core>
#include <Eigen/Geometry>
#include <numeric>
#include "sbs/physics/tetrahedral_mesh_boundary.h"
#include "sbs/topology/tetrahedron_set.h"

namespace sbs {
namespace physics {

/**
 * Surface mesh adapter
 */

tetrahedral_mesh_boundary_t::tetrahedral_mesh_boundary_t(topology::tetrahedron_set_t* mesh)
    : mesh_(mesh),
      vertex_index_map_{},
      tet_to_surface_vertex_index_map_{},
      triangle_index_map_{},
      vertices_{},
      triangles_{}
{
    extract_boundary_surface();
}

std::size_t tetrahedral_mesh_boundary_t::triangle_count() const
{
    return triangles_.size();
}

std::size_t tetrahedral_mesh_boundary_t::vertex_count() const
{
    return vertices_.size();
}

common::shared_vertex_surface_mesh_i::vertex_type
tetrahedral_mesh_boundary_t::vertex(std::size_t vi) const
{
    return vertices_[vi];
}

common::shared_vertex_surface_mesh_i::triangle_type
tetrahedral_mesh_boundary_t::triangle(std::size_t fi) const
{
    return triangles_[fi];
}

std::vector<index_type> const&
tetrahedral_mesh_boundary_t::surface_to_tetrahedral_mesh_index_map() const
{
    return vertex_index_map_;
}

index_type tetrahedral_mesh_boundary_t::from_surface_vertex(std::size_t vi) const
{
    return vertex_index_map_[vi];
}

index_type tetrahedral_mesh_boundary_t::from_surface_triangle(std::size_t fi) const
{
    return triangle_index_map_[fi];
}

void tetrahedral_mesh_boundary_t::extract_boundary_surface()
{
    vertex_index_map_.clear();
    triangle_index_map_.clear();
    tet_to_surface_vertex_index_map_.clear();
    vertices_.clear();
    triangles_.clear();

    std::size_t const tet_mesh_vertex_count = mesh_->vertex_count();
    std::size_t const tet_mesh_face_count   = mesh_->triangle_count();

    // maps tetrahedral mesh vertices to surface mesh vertices
    tet_to_surface_vertex_index_map_.resize(tet_mesh_vertex_count);

    // pre-allocate vertex storage heuristically, and triangle storage exactly
    vertex_index_map_.reserve(mesh_->triangle_count());
    triangle_index_map_.reserve(mesh_->triangle_count());
    triangles_.reserve(mesh_->triangle_count());
    vertices_.reserve(mesh_->vertex_count());

    std::vector<topology::triangle_t> const& triangles = mesh_->triangles();
    for (std::size_t fi = 0u; fi < triangles.size(); ++fi)
    {
        bool const is_interior_triangle = triangles[fi].incident_tetrahedron_indices().size() == 2u;
        if (is_interior_triangle)
            continue;

        triangle_index_map_.push_back(static_cast<index_type>(fi));

        topology::triangle_t const& boundary_triangle = triangles[fi];
        triangle_type surface_mesh_triangle{};

        for (std::size_t j = 0u; j < boundary_triangle.vertex_indices().size(); ++j)
        {
            index_type const vi = boundary_triangle.vertex_indices()[j];

            if (!tet_to_surface_vertex_index_map_[vi].has_value())
            {
                // update map
                index_type const new_vertex_index =
                    static_cast<index_type>(vertex_index_map_.size());
                tet_to_surface_vertex_index_map_[vi] = new_vertex_index;
                vertex_index_map_.push_back(vi);
                // Note: Vertex attributes should be set by the owning body in update_visual_model()
                vertices_.push_back({});
                surface_mesh_triangle.vertices[j] = new_vertex_index;
            }
            else
            {
                surface_mesh_triangle.vertices[j] = tet_to_surface_vertex_index_map_[vi].value();
            }
        }

        triangles_.push_back(surface_mesh_triangle);
    }
}

void tetrahedral_mesh_boundary_t::compute_normals()
{
    for (std::size_t i = 0u; i < triangles_.size(); ++i)
    {
        auto const v1 = triangles_[i].vertices[0u];
        auto const v2 = triangles_[i].vertices[1u];
        auto const v3 = triangles_[i].vertices[2u];

        Eigen::Vector3d const& p1 = vertices_[v1].position;
        Eigen::Vector3d const& p2 = vertices_[v2].position;
        Eigen::Vector3d const& p3 = vertices_[v3].position;

        Eigen::Vector3d const n = (p2 - p1).cross(p3 - p1);

        vertices_[v1].normal += n;
        vertices_[v2].normal += n;
        vertices_[v3].normal += n;
    }

    for (std::size_t i = 0u; i < vertices_.size(); ++i)
    {
        vertices_[i].normal.normalize();
    }
}

void tetrahedral_mesh_boundary_t::prepare_vertices_for_rendering()
{
    prepare_vertices_for_surface_rendering();
}

void tetrahedral_mesh_boundary_t::prepare_indices_for_rendering()
{
    prepare_indices_for_surface_rendering();
}

common::shared_vertex_surface_mesh_i::vertex_type&
tetrahedral_mesh_boundary_t::mutable_vertex(std::size_t vi)
{
    return vertices_[vi];
}

common::shared_vertex_surface_mesh_i::triangle_type&
tetrahedral_mesh_boundary_t::mutable_triangle(std::size_t f)
{
    return triangles_[f];
}

void tetrahedral_mesh_boundary_t::prepare_vertices_for_surface_rendering()
{
    std::size_t constexpr num_attributes_per_vertex = 9u;
    std::size_t const vertex_count                  = vertices_.size();
    std::vector<float> vertex_buffer{};
    vertex_buffer.reserve(vertex_count * num_attributes_per_vertex);

    for (vertex_type const& vertex : vertices_)
    {
        vertex_buffer.push_back(static_cast<float>(vertex.position.x()));
        vertex_buffer.push_back(static_cast<float>(vertex.position.y()));
        vertex_buffer.push_back(static_cast<float>(vertex.position.z()));

        // use triangle normal since face-based
        vertex_buffer.push_back(static_cast<float>(vertex.normal.x()));
        vertex_buffer.push_back(static_cast<float>(vertex.normal.y()));
        vertex_buffer.push_back(static_cast<float>(vertex.normal.z()));

        vertex_buffer.push_back(vertex.color.x());
        vertex_buffer.push_back(vertex.color.y());
        vertex_buffer.push_back(vertex.color.z());
    }

    transfer_vertices_for_rendering(std::move(vertex_buffer));
}

void tetrahedral_mesh_boundary_t::prepare_indices_for_surface_rendering()
{
    std::size_t const index_count = triangles_.size() * 3u;
    std::vector<std::uint32_t> index_buffer{};
    index_buffer.reserve(index_count);
    for (triangle_type const& f : triangles_)
    {
        index_buffer.push_back(f.vertices[0]);
        index_buffer.push_back(f.vertices[1]);
        index_buffer.push_back(f.vertices[2]);
    }

    transfer_indices_for_rendering(std::move(index_buffer));
}

void tetrahedral_mesh_boundary_t::prepare_vertices_for_wireframe_rendering()
{
    // no-op
}

void tetrahedral_mesh_boundary_t::prepare_indices_for_wireframe_rendering()
{
    // no-op
}

topology::tetrahedron_set_t const* tetrahedral_mesh_boundary_t::tetrahedral_mesh() const
{
    return mesh_;
}

topology::tetrahedron_set_t* tetrahedral_mesh_boundary_t::tetrahedral_mesh()
{
    return mesh_;
}

} // namespace physics
} // namespace sbs