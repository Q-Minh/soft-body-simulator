#ifndef SBS_PHYSICS_MESH_H
#define SBS_PHYSICS_MESH_H

#include "common/mesh.h"

#include <Eigen/Core>
#include <array>
#include <memory>
#include <optional>
#include <vector>

namespace sbs {
namespace physics {

using index_type  = std::uint32_t;
using scalar_type = double;

class tetrahedron_t;
class edge_t;
class triangle_t;

class vertex_t
{
  public:
    vertex_t();
    vertex_t(vertex_t const&) = default;
    vertex_t(vertex_t&&)      = default;

    using position_type = Eigen::Vector3d;
    using velocity_type = Eigen::Vector3d;
    using force_type    = Eigen::Vector3d;
    using normal_type   = Eigen::Vector3d;
    using color_type    = Eigen::Vector3f;

    position_type const& position() const;
    position_type& position();

    velocity_type const& velocity() const;
    velocity_type& velocity();

    force_type const& force() const;
    force_type& force();

    scalar_type const& mass() const;
    scalar_type& mass();

    normal_type const& normal() const;
    normal_type& normal();

    color_type const& color() const;
    color_type& color();

    bool fixed() const;
    bool& fixed();

  private:
    position_type p_;
    velocity_type v_;
    force_type f_;
    normal_type n_;
    double m_;
    bool fixed_;
    color_type c_;
};

class edge_t
{
  public:
    edge_t() = default;
    edge_t(index_type v1, index_type v2);
    edge_t(edge_t const&) = default;
    edge_t(edge_t&&)      = default;

    index_type const& v1() const;
    index_type const& v2() const;

    index_type& v1();
    index_type& v2();

    std::vector<index_type>& adjacent_tetrahedron_indices();
    std::vector<index_type> const& adjacent_tetrahedron_indices() const;

    std::vector<index_type>& adjacent_triangle_indices();
    std::vector<index_type> const& adjacent_triangle_indices() const;

    /**
     * @brief Checks if other edge is the same as this edge, disregarding orientation.
     * @param other Edge to test against
     * @return True if other edge has same vertices as this edge
     */
    bool operator==(edge_t const& other) const;

    /**
     * @brief
     * Specifies ordering of edge in a collection of edges for use in a std::map, for
     * example.
     * @param other Edge to compare against.
     * @return True if this edge comes before other edge in the ordering.
     */
    bool operator<(edge_t const& other) const;

    /**
     * @brief
     * Checks if other edge has opposite orientation from this edge.
     * @pre (*this == other) is true.
     * @param other The edge to test against
     * @return True if other has opposite orientation from this edge
     */
    bool is_reverse_of(edge_t const& other) const;

  private:
    std::array<index_type, 2u> v_;
    std::vector<index_type> adjacent_tets_;
    std::vector<index_type> adjacent_triangles_;
};

class triangle_t
{
  public:
    triangle_t() = default;
    triangle_t(index_type v1, index_type v2, index_type v3);
    triangle_t(triangle_t const&);
    triangle_t(triangle_t&&) = default;

    index_type const& v1() const;
    index_type const& v2() const;
    index_type const& v3() const;

    index_type& v1();
    index_type& v2();
    index_type& v3();

    /**
     * @brief
     * Return list of edge copies of this triangle's edges containing vertex indices, but no
     * connectivity information. Can be used even if adjacency information for this triangle
     * has not been built.
     * @return Edge copies of this triangle
     */
    std::array<edge_t, 3u> edges_copy() const;

    /**
     * @brief Checks if other triangle is the same as this triangle, disregarding orientation.
     * @param other Triangle to test against
     * @return True if other triangle has same vertices as this triangle
     */
    bool operator==(triangle_t const& other) const;

    /**
     * @brief
     * Specifies ordering of triangle in a collection of triangles for use in a std::map, for
     * example.
     * @param other Triangle to compare against.
     * @return True if this triangle comes before other triangle in the ordering.
     */
    bool operator<(triangle_t const& other) const;

    /**
     * @brief
     * Checks if other triangle has opposite orientation from this triangle.
     * @pre (*this == other) is true.
     * @param other The triangle to test against
     * @return True if other has opposite orientation from this triangle
     */
    bool is_reverse_of(triangle_t const& other) const;

    bool is_boundary_triangle() const;
    bool is_interior_triangle() const;

    std::vector<index_type> const& adjacent_tetrahedron_indices() const;
    std::vector<index_type>& adjacent_tetrahedron_indices();

    std::unique_ptr<std::array<index_type, 3u>> const& adjacent_edge_indices() const;
    std::unique_ptr<std::array<index_type, 3u>>& adjacent_edge_indices();

  private:
    std::array<index_type, 3u> v_;
    std::vector<index_type> adjacent_tets_;
    std::unique_ptr<std::array<index_type, 3u>> edges_;
};

class tetrahedron_t
{
  public:
    tetrahedron_t() = default;
    tetrahedron_t(index_type v1, index_type v2, index_type v3, index_type v4);
    tetrahedron_t(tetrahedron_t const& other);
    tetrahedron_t(tetrahedron_t&&) = default;

    index_type const& v1() const;
    index_type const& v2() const;
    index_type const& v3() const;
    index_type const& v4() const;

    index_type& v1();
    index_type& v2();
    index_type& v3();
    index_type& v4();

    scalar_type const& mass_density() const;
    scalar_type& mass_density();

    /**
     * @brief
     * Return list of edge copies of this tetrahedron's edges containing vertex indices, but no
     * connectivity information. Can be used even if adjacency information for this tetrahedron
     * has not been built.
     * @return Edge copies of this tetrahedron
     */
    std::array<edge_t, 6u> edges_copy() const;

    /**
     * @brief
     * Return list of face copies of this tetrahedron's faces containing vertex indices, but no
     * connectivity information. Can be used even if adjacency information for this tetrahedron
     * has not been built.
     * @return Face copies of this tetrahedron
     */
    std::array<triangle_t, 4u> faces_copy() const;

    std::unique_ptr<std::array<index_type, 6u>> const& edge_indices() const;
    std::unique_ptr<std::array<index_type, 6u>>& edge_indices();

    std::unique_ptr<std::array<index_type, 4u>> const& face_indices() const;
    std::unique_ptr<std::array<index_type, 4u>>& face_indices();

  private:
    std::array<index_type, 4u> v_;
    std::unique_ptr<std::array<index_type, 6u>> edges_;
    std::unique_ptr<std::array<index_type, 4u>> faces_;
    double rho_; ///< mass density
};

struct build_topology_parameters_t
{
    bool vertex_to_edge{
        false}; ///< Vertices must have references to their adjacent edges. Not supported yet.
    bool vertex_to_triangle{
        false}; ///< Vertices must have references to their adjacent triangles. Not supported yet.
    bool vertex_to_tetrahedra{
        false}; ///< Vertices must have references to their adjacent tetrahedra. Not supported yet.

    bool triangle_to_edge{false}; ///< Triangles must have references to their edges
    bool triangle_to_tetrahedra{
        false}; ///< Triangles must have references to their adjacent tetrahedra

    bool tetrahedron_to_edge{false};     ///< Tetrahedra must have references to their edges
    bool tetrahedron_to_triangle{false}; ///< Tetrahedra must have references to their faces

    bool edge_to_triangle{false};   ///< Edges must have references to their adjacent triangles
    bool edge_to_tetrahedra{false}; ///< Edges must have references to their adjacent tetrahedra

    void include_vertex_to_edge_adjacency();
    void include_vertex_to_triangle_adjacency();
    void include_vertex_to_tetrahedra_adjacency();
    void include_triangle_to_edge_adjacency();
    void include_triangle_to_tetrahedra_adjacency();
    void include_tetrahedron_to_edge_adjacency();
    void include_tetrahedron_to_triangle_adjacency();
    void include_edge_to_triangle_adjacency();
    void include_edge_to_tetrahedra_adjacency();

    void exclude_vertex_to_edge_adjacency();
    void exclude_vertex_to_triangle_adjacency();
    void exclude_vertex_to_tetrahedra_adjacency();
    void exclude_triangle_to_edge_adjacency();
    void exclude_triangle_to_tetrahedra_adjacency();
    void exclude_tetrahedron_to_edge_adjacency();
    void exclude_tetrahedron_to_triangle_adjacency();
    void exclude_edge_to_triangle_adjacency();
    void exclude_edge_to_tetrahedra_adjacency();
};

class topological_simulated_tetrahedral_mesh_t
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
    void build_topology(build_topology_parameters_t const& params);

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
    renderable_topological_simulated_tetrahedral_mesh_t(
        common::geometry_t const& geometry,
        build_topology_parameters_t const& params = build_topology_parameters_t{});

    // TODO: Support wireframe rendering
    virtual void prepare_vertices_for_rendering() override;
    virtual void prepare_indices_for_rendering() override;
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

#endif // SBS_PHYSICS_MESH_H