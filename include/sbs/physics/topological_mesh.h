#ifndef SBS_PHYSICS_TOPOLOGICAL_MESH_H
#define SBS_PHYSICS_TOPOLOGICAL_MESH_H

#include "sbs/physics/mesh.h"

namespace sbs {
namespace physics {

class topological_simulated_tetrahedral_mesh_t
    : public simulated_mesh_i<topological_simulated_tetrahedral_mesh_t>
{
  public:
    topological_simulated_tetrahedral_mesh_t() = default;

    /**
     * @brief
     * Builds the topological mesh according to the given geometry
     * and topological information requested. The built vertices will
     * take the geometry's positions, normals and colors. The built
     * tetrahedra will take the geometry's indices. Topological
     * information is then built according to params.
     * @param geometry The geometry of the mesh
     * @param params The required topological information to build on this mesh
     */
    topological_simulated_tetrahedral_mesh_t(
        common::geometry_t const& geometry,
        build_topology_parameters_t const& params = build_topology_parameters_t{});

    std::vector<vertex_t> const& vertices() const;
    std::vector<vertex_t>& vertices();

    std::vector<edge_t> const& edges() const;
    std::vector<edge_t>& edges();

    std::vector<triangle_t> const& triangles() const;
    std::vector<triangle_t>& triangles();

    std::vector<tetrahedron_t> const& tetrahedra() const;
    std::vector<tetrahedron_t>& tetrahedra();

    build_topology_parameters_t const& topology_parameters() const;
    build_topology_parameters_t& topology_parameters();
    void build_topology();

    std::vector<triangle_t const*> boundary_triangles() const;

    std::size_t boundary_triangle_count() const;

    /**
     * @brief Get adjacent tetrahedra to the given edge
     * @param edge
     * @return
     */
    std::vector<tetrahedron_t const*> edge_adjacent_tetrahedra(index_type edge) const;

    /**
     * @brief Get adjacent tetrahedra to the given triangle
     * @param triangle
     * @return
     */
    std::vector<tetrahedron_t const*> triangle_adjacent_tetrahedra(index_type triangle) const;

    /**
     * @brief Get edges of the given tetrahedron
     * @param tetrahedron
     * @return
     */
    std::array<edge_t const*, 6u> edges_of_tetrahedron(index_type tetrahedron) const;

    /**
     * @brief Get faces of the given tetrahedron
     * @param tetrahedron
     * @return
     */
    std::array<triangle_t const*, 4u> faces_of_tetrahedron(index_type tetrahedron) const;

    /**
     * @brief Get vertices of the given tetrahedron
     * @param tetrahedron
     * @return
     */
    std::array<vertex_t const*, 4u> vertices_of_tetrahedron(index_type tetrahedron) const;

    /**
     * @brief Get adjacent tetrahedra to the given edge
     * @param edge
     * @return
     */
    std::vector<tetrahedron_t*> edge_adjacent_tetrahedra(index_type edge);

    /**
     * @brief Get adjacent tetrahedra to the given triangle
     * @param triangle
     * @return
     */
    std::vector<tetrahedron_t*> triangle_adjacent_tetrahedra(index_type triangle);

    /**
     * @brief Get edges of the given tetrahedron
     * @param tetrahedron
     * @return
     */
    std::array<edge_t*, 6u> edges_of_tetrahedron(index_type tetrahedron);

    /**
     * @brief Get faces of the given tetrahedron
     * @param tetrahedron
     * @return
     */
    std::array<triangle_t*, 4u> faces_of_tetrahedron(index_type tetrahedron);

    /**
     * @brief Get vertices of the given tetrahedron
     * @param tetrahedron
     * @return
     */
    std::array<vertex_t*, 4u> vertices_of_tetrahedron(index_type tetrahedron);

    bool is_boundary_edge(index_type edge) const;
    bool is_interior_edge(index_type edge) const;
    bool is_boundary_triangle(index_type triangle) const;
    bool is_interior_triangle(index_type triangle) const;

  private:
    std::vector<vertex_t> vertices_;
    std::vector<edge_t> edges_;
    std::vector<triangle_t> triangles_;
    std::vector<tetrahedron_t> tetrahedra_;

    build_topology_parameters_t topology_params_;
};

/**
 * @brief
 * Mesh that can be used for simulation as well as rendering.
 * Renders as a face-based mesh.
 */
class renderable_topological_simulated_tetrahedral_mesh_t
    : public topological_simulated_tetrahedral_mesh_t,
      public common::renderable_node_t
{
  public:
    renderable_topological_simulated_tetrahedral_mesh_t() = default;

    renderable_topological_simulated_tetrahedral_mesh_t(
        common::geometry_t const& geometry,
        build_topology_parameters_t const& params = build_topology_parameters_t{});

    virtual void prepare_vertices_for_rendering() override;
    virtual void prepare_indices_for_rendering() override;

  private:
    void prepare_vertices_for_surface_rendering();
    void prepare_indices_for_surface_rendering();
    void prepare_vertices_for_wireframe_rendering();
    void prepare_indices_for_wireframe_rendering();
};

/**
 * @brief Surface mesh representation of a tetrahedral mesh's boundary surface.
 */
class tetrahedral_mesh_surface_mesh_adapter_t : public common::shared_vertex_surface_mesh_i
{
  public:
    using vertex_type   = common::shared_vertex_surface_mesh_i::vertex_type;
    using triangle_type = common::shared_vertex_surface_mesh_i::triangle_type;

    tetrahedral_mesh_surface_mesh_adapter_t(topological_simulated_tetrahedral_mesh_t* mesh);

    virtual std::size_t triangle_count() const override;
    virtual std::size_t vertex_count() const override;

    virtual shared_vertex_surface_mesh_i::vertex_type vertex(std::size_t vi) const override;
    virtual shared_vertex_surface_mesh_i::triangle_type triangle(std::size_t f) const override;

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

    topological_simulated_tetrahedral_mesh_t const* tetrahedral_mesh() const;
    topological_simulated_tetrahedral_mesh_t* tetrahedral_mesh();

  private:
    topological_simulated_tetrahedral_mesh_t* mesh_;
    std::vector<index_type>
        vertex_index_map_; ///< Maps from surface vertex indices to tet vertex indices
    std::vector<std::optional<index_type>>
        tet_to_surface_vertex_index_map_; ///< Maps from tet vertex indices to surface vertex
                                          ///< indices
    std::vector<index_type>
        triangle_index_map_; ///< Maps from surface face indices to tet triangle indices
    std::vector<Eigen::Vector3d> vertex_normals_; ///< Vertex normals of boundary vertices
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_TOPOLOGICAL_MESH_H