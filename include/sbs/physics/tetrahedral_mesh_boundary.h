#ifndef SBS_PHYSICS_TETRAHEDRAL_MESH_BOUNDARY_H
#define SBS_PHYSICS_TETRAHEDRAL_MESH_BOUNDARY_H

#include <optional>
#include <sbs/aliases.h>
#include <sbs/common/mesh.h>
#include <sbs/common/node.h>
#include <vector>

namespace sbs {
namespace physics {

// Forward declares
class tetrahedron_set_t;

/**
 * @brief Surface mesh representation of a tetrahedral mesh's boundary surface.
 */
class tetrahedral_mesh_boundary_t : public common::shared_vertex_surface_mesh_i
{
  public:
    using vertex_type   = common::shared_vertex_surface_mesh_i::vertex_type;
    using triangle_type = common::shared_vertex_surface_mesh_i::triangle_type;

    tetrahedral_mesh_boundary_t() = default;
    tetrahedral_mesh_boundary_t(tetrahedron_set_t* mesh);

    tetrahedral_mesh_boundary_t(tetrahedral_mesh_boundary_t const& other) = default;
    tetrahedral_mesh_boundary_t(tetrahedral_mesh_boundary_t&& other)      = default;
    tetrahedral_mesh_boundary_t& operator=(tetrahedral_mesh_boundary_t const& other) = default;
    tetrahedral_mesh_boundary_t& operator=(tetrahedral_mesh_boundary_t&& other) = default;

    virtual std::size_t triangle_count() const override;
    virtual std::size_t vertex_count() const override;

    virtual shared_vertex_surface_mesh_i::vertex_type vertex(std::size_t vi) const override;
    virtual shared_vertex_surface_mesh_i::triangle_type triangle(std::size_t f) const override;

    virtual void prepare_vertices_for_rendering() override;
    virtual void prepare_indices_for_rendering() override;

    /**
     * @brief Gets the mapping from this surface mesh's vertex indices to its underlying tetrahedral
     * mesh's vertex indices
     * @return The map from surface mesh vertex indices to tetrahedral mesh indices
     */
    std::vector<index_type> const& surface_to_tetrahedral_mesh_index_map() const;

    /**
     * @brief Returns the tetrahedral mesh's corresponding index of the surface mesh's vertex index
     * vi.
     * @param vi The surface mesh's vertex index
     * @return The tetrahedral mesh's vertex index corresponding to the surface mesh's vertex index
     * vi
     */
    index_type from_surface_vertex(std::size_t vi) const;

    /**
     * @brief Returns the tetrahedral mesh's corresponding index of the surface mesh's triangle
     * index fi.
     * @param fi The surface mesh's triangle index
     * @return The tetrahedral mesh's triangle index corresponding to the surface mesh's triangle
     * index
     */
    index_type from_surface_triangle(std::size_t fi) const;

    /**
     * @brief Recompute the tetrahedral mesh's boundary surface mesh
     */
    void extract_boundary_surface();

    void compute_normals();

    tetrahedron_set_t const* tetrahedral_mesh() const;
    tetrahedron_set_t* tetrahedral_mesh();

  private:
    void prepare_vertices_for_surface_rendering();
    void prepare_indices_for_surface_rendering();
    void prepare_vertices_for_wireframe_rendering();
    void prepare_indices_for_wireframe_rendering();

    tetrahedron_set_t* mesh_;
    std::vector<index_type>
        vertex_index_map_; ///< Maps from surface vertex indices to tet vertex indices
    std::vector<std::optional<index_type>>
        tet_to_surface_vertex_index_map_; ///< Maps from tet vertex indices to surface vertex
                                          ///< indices
    std::vector<index_type>
        triangle_index_map_; ///< Maps from surface face indices to tet triangle indices

    std::vector<vertex_type> vertices_;    ///< Surface mesh vertices
    std::vector<triangle_type> triangles_; ///< Surface mesh triangles
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_TETRAHEDRAL_MESH_BOUNDARY_H