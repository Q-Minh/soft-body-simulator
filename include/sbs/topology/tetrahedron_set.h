#ifndef SBS_TOPOLOGY_TETRAHEDRON_SET_H
#define SBS_TOPOLOGY_TETRAHEDRON_SET_H

#include <array>
#include <map>
#include <sbs/aliases.h>
#include <set>
#include <vector>

namespace sbs {
namespace topology {

class vertex_t
{
  public:
    vertex_t() = default;

    vertex_t(index_type vi);

    index_type vi() const;

    std::vector<index_type> const& incident_edge_indices() const;
    std::vector<index_type> const& incident_triangle_indices() const;
    std::vector<index_type> const& incident_tetrahedron_indices() const;
    std::vector<index_type>& incident_edge_indices();
    std::vector<index_type>& incident_triangle_indices();
    std::vector<index_type>& incident_tetrahedron_indices();
    bool operator==(vertex_t const& other) const;
    bool operator!=(vertex_t const& other) const;
    bool operator<(vertex_t const& other) const;

  private:
    index_type vi_;                              ///< Index of this vertex
    std::vector<index_type> incident_edges_;     ///< Indices to incident edges
    std::vector<index_type> incident_triangles_; ///< Indices to incident triangles
    std::vector<index_type> incident_tets_;      ///< Indices to incident tetrahedra
};

class edge_t
{
  public:
    edge_t() = default;
    edge_t(index_type v1, index_type v2);
    edge_t(vertex_t const& v1, vertex_t const& v2);
    std::array<index_type, 2u> const& vertex_indices() const;
    index_type v1() const;
    index_type v2() const;
    std::vector<index_type> const& incident_triangle_indices() const;
    std::vector<index_type> const& incident_tetrahedron_indices() const;
    std::vector<index_type>& incident_triangle_indices();
    std::vector<index_type>& incident_tetrahedron_indices();
    std::uint8_t id_of_vertex(index_type vi) const;
    bool is_reverse_of(edge_t const& other) const;
    bool operator==(edge_t const&) const;
    bool operator<(edge_t const&) const;

  private:
    std::array<index_type, 2u> v_;               ///< Vertex indices
    std::vector<index_type> incident_triangles_; ///< Indices to incident triangles
    std::vector<index_type> incident_tets_;      ///< indices to incident tetrahedra
};

class triangle_t
{
  public:
    triangle_t() = default;
    triangle_t(index_type v1, index_type v2, index_type v3);
    triangle_t(vertex_t const& v1, vertex_t const& v2, vertex_t const& v3);
    std::array<index_type, 3u> const& vertex_indices() const;
    index_type v1() const;
    index_type v2() const;
    index_type v3() const;

    std::array<edge_t, 3u> edges_copy() const;
    std::array<index_type, 3u> const& edge_indices() const;
    std::array<index_type, 3u>& edge_indices();

    std::vector<index_type> const& incident_tetrahedron_indices() const;
    std::vector<index_type>& incident_tetrahedron_indices();

    std::uint8_t id_of_vertex(index_type vi) const;
    std::uint8_t id_of_edge(index_type ei) const;
    bool is_reverse_of(triangle_t const& other) const;

    bool operator==(triangle_t const& other) const;
    bool operator<(triangle_t const& other) const;

  private:
    std::array<index_type, 3u> v_;          ///< Vertex indices
    std::array<index_type, 3u> edges_;      ///< Indices to this triangle's edges
    std::vector<index_type> incident_tets_; ///< Indices to tetahedra incident on this triangle
};

class tetrahedron_t
{
  public:
    tetrahedron_t() = default;
    tetrahedron_t(index_type v1, index_type v2, index_type v3, index_type v4);
    tetrahedron_t(vertex_t const& v1, vertex_t const& v2, vertex_t const& v3, vertex_t const& v4);

    std::array<index_type, 4u> const& vertex_indices() const;
    std::array<index_type, 4u>& vertex_indices();
    index_type v1() const;
    index_type v2() const;
    index_type v3() const;
    index_type v4() const;

    index_type& v1();
    index_type& v2();
    index_type& v3();
    index_type& v4();
    std::array<edge_t, 6u> edges_copy() const;
    std::array<triangle_t, 4u> faces_copy() const;
    std::array<index_type, 6u> const& edge_indices() const;
    std::array<index_type, 4u> const& face_indices() const;
    std::array<index_type, 6u>& edge_indices();
    std::array<index_type, 4u>& face_indices();
    std::uint8_t id_of_vertex(index_type vi) const;
    std::uint8_t id_of_edge(index_type ei) const;
    std::uint8_t id_of_face(index_type fi) const;
    bool operator==(tetrahedron_t const& other) const;
    bool operator<(tetrahedron_t const& other) const;

  private:
    std::array<index_type, 4u> v_{};   ///< Vertex indices
    std::array<index_type, 6u> edges_; ///< Edge indices
    std::array<index_type, 4u> faces_; ///< Face indices
};

class vertex_set_t
{
  public:
    vertex_set_t() = default;
    template <class VertexIterator>
    vertex_set_t(VertexIterator begin, VertexIterator end);
    vertex_t const& vertex(std::size_t vi) const;
    vertex_t& vertex(std::size_t vi);

    bool contains_vertex(vertex_t const& v) const;

    index_type add_vertex();
    void add_vertex(vertex_t const& vertex);
    void add_vertex(index_type const vi);
    std::size_t vertex_count() const;
    void reserve_vertices(std::size_t count);
    void clear();
    std::vector<vertex_t> const& vertices() const;

    void remove_vertex_to_edge_incidency(index_type const vi, index_type ei);
    void remove_vertex_to_triangle_incidency(index_type const vi, index_type fi);
    void remove_vertex_to_tetrahedron_incidency(index_type const vi, index_type ti);
    bool operator==(vertex_set_t const& other) const;

  private:
    std::vector<vertex_t> vertices_; ///< Vertices of this set
};

template <class VertexIterator>
inline vertex_set_t::vertex_set_t(VertexIterator begin, VertexIterator end) : vertices_{begin, end}
{
}

class edge_set_t : public vertex_set_t
{
  public:
    edge_set_t() = default;

    edge_t const& edge(index_type ei) const;
    edge_t& edge(index_type ei);
    index_type ei(edge_t const& edge) const;
    index_type add_edge(edge_t const& edge);
    edge_t remove_edge(index_type ei);
    edge_t remove_edge(edge_t const& edge);
    std::size_t edge_count() const;
    bool contains_edge(edge_t const& edge) const;
    void reserve_edges(std::size_t count);
    void clear();
    bool is_safe_to_iterate_over_edges() const;

    std::vector<edge_t> const& edges() const;
    std::vector<edge_t>& edges();

    std::map<edge_t, index_type>::const_iterator safe_edges_begin() const;
    std::map<edge_t, index_type>::const_iterator safe_edges_end() const;

    void remove_edge_to_triangle_incidency(index_type const ei, index_type fi);
    void remove_edge_to_tetrahedron_incidency(index_type const ei, index_type ti);
    void create_vertex_to_edge_incidency(index_type const ei);
    bool operator==(edge_set_t const& other) const;

  protected:
    std::map<edge_t, index_type> edge_map_; ///< Map from edge to edge indices
    std::set<index_type>
        edge_garbage_collector_; ///< Set of indices that can be reused for edge creation

  private:
    std::vector<edge_t> edges_; ///< Edges of this set
};

class triangle_set_t : public edge_set_t
{
  public:
    triangle_set_t() = default;

    triangle_t const& triangle(index_type fi) const;
    triangle_t& triangle(index_type fi);
    index_type fi(triangle_t const& triangle) const;

    void add_triangle(triangle_t const& triangle);
    triangle_t remove_triangle(index_type fi);
    triangle_t remove_triangle(triangle_t const& triangle);

    std::size_t triangle_count() const;
    bool contains_triangle(triangle_t const& triangle) const;
    void reserve_triangles(std::size_t count);
    void clear();
    bool is_safe_to_iterate_over_triangles() const;

    std::vector<triangle_t> const& triangles() const;
    std::vector<triangle_t>& triangles();
    std::map<triangle_t, index_type>::const_iterator safe_triangles_begin() const;
    std::map<triangle_t, index_type>::const_iterator safe_triangles_end() const;

    void create_vertex_to_triangle_incidency(index_type const fi);
    void create_edge_to_triangle_incidency(index_type const fi);
    void remove_triangle_to_tetrahedron_incidency(index_type const fi, index_type const ti);
    bool operator==(triangle_set_t const& other) const;

  protected:
    std::map<triangle_t, index_type> triangle_map_;   ///< Map of triangles to their indices
    std::set<index_type> triangle_garbage_collector_; ///< Indices of deleted triangles

  private:
    std::vector<triangle_t> triangles_; ///< Triangles of this set
};

class tetrahedron_set_t : public triangle_set_t
{
  public:
    tetrahedron_set_t() = default;

    tetrahedron_t const& tetrahedron(index_type ti) const;
    tetrahedron_t& tetrahedron(index_type ti);

    void add_tetrahedron(tetrahedron_t const& tetrahedron);
    tetrahedron_t remove_tetrahedron(index_type ti);

    std::size_t tetrahedron_count() const;
    void reserve_tetrahedra(std::size_t count);
    void clear();
    bool is_safe_to_iterate_over_tetrahedra() const;

    std::vector<tetrahedron_t> const& tetrahedra() const;
    std::vector<tetrahedron_t>& tetrahedra();

    void create_vertex_to_tetrahedron_incidency(index_type const ti);
    void create_edge_to_tetrahedron_incidency(index_type const ti);
    void create_triangle_to_tetrahedron_incidency(index_type const ti);

    void collect_garbage();

    bool is_boundary_tetrahedron(index_type const ti) const;
    bool is_boundary_triangle(index_type const fi) const;
    bool is_boundary_edge(index_type const ei) const;
    bool is_boundary_vertex(index_type const vi) const;

    bool is_boundary_tetrahedron(tetrahedron_t const& t) const;
    bool is_boundary_triangle(triangle_t const& f) const;
    bool is_boundary_edge(edge_t const& e) const;
    bool is_boundary_vertex(vertex_t const& v) const;

    std::vector<tetrahedron_t> boundary_tetrahedra() const;
    std::vector<index_type> boundary_tetrahedron_indices() const;
    std::vector<triangle_t> boundary_triangles() const;
    std::vector<index_type> boundary_triangle_indices() const;
    std::vector<edge_t> boundary_edges() const;
    std::vector<index_type> boundary_edge_indices() const;
    std::vector<vertex_t> boundary_vertices() const;
    std::vector<index_type> boundary_vertex_indices() const;

    std::vector<tetrahedron_t> interior_tetrahedra() const;
    std::vector<index_type> interior_tetrahedron_indices() const;
    std::vector<triangle_t> interior_triangles() const;
    std::vector<index_type> interior_triangle_indices() const;
    std::vector<edge_t> interior_edges() const;
    std::vector<index_type> interior_edge_indices() const;
    std::vector<vertex_t> interior_vertices() const;
    std::vector<index_type> interior_vertex_indices() const;

    std::vector<triangle_t> oriented_boundary_triangles() const;

    bool operator==(tetrahedron_set_t const& other) const;

  protected:
    std::set<index_type> tetrahedron_garbage_collector_; ///< Indices of deleted tetrahedra

  private:
    void swap_edges(index_type ei, index_type eip);
    void swap_triangles(index_type fi, index_type fip);
    void swap_tetrahedra(index_type ti, index_type tip);

    std::vector<tetrahedron_t> tetrahedra_; ///< Tetrahedra of this set
};

} // namespace topology
} // namespace sbs

#endif // SBS_TOPOLOGY_TETRAHEDRON_SET_H