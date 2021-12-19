#include <algorithm>
#include <cassert>
#include <iterator>
#include <sbs/topology/tetrahedron_set.h>

namespace sbs {
namespace topology {

/**
 * Vertex implementation
 */

vertex_t::vertex_t(index_type vi)
    : vi_(vi), incident_edges_{}, incident_triangles_{}, incident_tets_{}
{
}

index_type vertex_t::vi() const
{
    return vi_;
}

std::vector<index_type> const& vertex_t::incident_edge_indices() const
{
    return incident_edges_;
}

std::vector<index_type> const& vertex_t::incident_triangle_indices() const
{
    return incident_triangles_;
}

std::vector<index_type> const& vertex_t::incident_tetrahedron_indices() const
{
    return incident_tets_;
}

std::vector<index_type>& vertex_t::incident_edge_indices()
{
    return incident_edges_;
}

std::vector<index_type>& vertex_t::incident_triangle_indices()
{
    return incident_triangles_;
}

std::vector<index_type>& vertex_t::incident_tetrahedron_indices()
{
    return incident_tets_;
}

bool vertex_t::operator==(vertex_t const& other) const
{
    return vi_ == other.vi_;
}

bool vertex_t::operator!=(vertex_t const& other) const
{
    return vi_ != other.vi_;
}

bool vertex_t::operator<(vertex_t const& other) const
{
    return vi_ < other.vi_;
}

/**
 * Edge implementation
 */

edge_t::edge_t(index_type v1, index_type v2) : v_{v1, v2}, incident_triangles_{}, incident_tets_{}
{
}

edge_t::edge_t(vertex_t const& v1, vertex_t const& v2)
    : v_{v1.vi(), v2.vi()}, incident_triangles_{}, incident_tets_{}
{
}

std::array<index_type, 2u> const& edge_t::vertex_indices() const
{
    return v_;
}

index_type edge_t::v1() const
{
    return v_[0];
}

index_type edge_t::v2() const
{
    return v_[1];
}

std::vector<index_type> const& edge_t::incident_triangle_indices() const
{
    return incident_triangles_;
}

std::vector<index_type> const& edge_t::incident_tetrahedron_indices() const
{
    return incident_tets_;
}

std::vector<index_type>& edge_t::incident_triangle_indices()
{
    return incident_triangles_;
}

std::vector<index_type>& edge_t::incident_tetrahedron_indices()
{
    return incident_tets_;
}

std::uint8_t edge_t::id_of_vertex(index_type vi) const
{
    return static_cast<std::uint8_t>(
        std::distance(v_.begin(), std::find(v_.begin(), v_.end(), vi)));
}

bool edge_t::is_reverse_of(edge_t const& other) const
{
    return v1() == other.v2() && v2() == other.v1();
}

bool edge_t::operator==(edge_t const& other) const
{
    std::array<index_type, 2u> v      = v_;
    std::array<index_type, 2u> vother = other.v_;

    std::sort(v.begin(), v.end());
    std::sort(vother.begin(), vother.end());

    return v == vother;
}

bool edge_t::operator<(edge_t const& other) const
{
    std::array<index_type, 2u> v      = v_;
    std::array<index_type, 2u> vother = other.v_;

    std::sort(v.begin(), v.end());
    std::sort(vother.begin(), vother.end());

    return v < vother;
}

/**
 * Triangle implementation
 */

triangle_t::triangle_t(index_type v1, index_type v2, index_type v3)
    : v_{v1, v2, v3}, edges_{}, incident_tets_{}
{
}

triangle_t::triangle_t(vertex_t const& v1, vertex_t const& v2, vertex_t const& v3)
    : v_{v1.vi(), v2.vi(), v3.vi()}, edges_{}, incident_tets_{}
{
}

std::array<index_type, 3u> const& triangle_t::vertex_indices() const
{
    return v_;
}

index_type triangle_t::v1() const
{
    return v_[0];
}

index_type triangle_t::v2() const
{
    return v_[1];
}

index_type triangle_t::v3() const
{
    return v_[2];
}

std::array<edge_t, 3u> triangle_t::edges_copy() const
{
    return std::array<edge_t, 3u>{edge_t{v1(), v2()}, edge_t{v2(), v3()}, edge_t{v3(), v1()}};
}

std::array<index_type, 3u> const& triangle_t::edge_indices() const
{
    return edges_;
}

std::vector<index_type> const& triangle_t::incident_tetrahedron_indices() const
{
    return incident_tets_;
}

std::array<index_type, 3u>& triangle_t::edge_indices()
{
    return edges_;
}

std::vector<index_type>& triangle_t::incident_tetrahedron_indices()
{
    return incident_tets_;
}

std::uint8_t triangle_t::id_of_vertex(index_type vi) const
{
    return static_cast<std::uint8_t>(
        std::distance(v_.begin(), std::find(v_.begin(), v_.end(), vi)));
}

std::uint8_t triangle_t::id_of_edge(index_type ei) const
{
    return static_cast<std::uint8_t>(
        std::distance(edges_.begin(), std::find(edges_.begin(), edges_.end(), ei)));
}

bool triangle_t::is_reverse_of(triangle_t const& other) const
{
    std::array<index_type, 3u> v      = v_;
    std::array<index_type, 3u> vother = other.v_;

    auto const v_min      = std::min_element(v.begin(), v.end());
    auto const vother_min = std::min_element(vother.begin(), vother.end());

    std::rotate(v.begin(), v_min, v.end());
    std::rotate(vother.begin(), vother_min, vother.end());

    return v != vother;
}

bool triangle_t::operator==(triangle_t const& other) const
{
    std::array<index_type, 3u> v      = v_;
    std::array<index_type, 3u> vother = other.v_;

    std::sort(v.begin(), v.end());
    std::sort(vother.begin(), vother.end());

    return v == vother;
}

bool triangle_t::operator<(triangle_t const& other) const
{
    std::array<index_type, 3u> v      = v_;
    std::array<index_type, 3u> vother = other.v_;

    std::sort(v.begin(), v.end());
    std::sort(vother.begin(), vother.end());

    return v < vother;
}

/**
 * Tetrahedron implementation
 */

tetrahedron_t::tetrahedron_t(index_type v1, index_type v2, index_type v3, index_type v4)
    : v_{v1, v2, v3, v4}, edges_{}, faces_{}
{
}

tetrahedron_t::tetrahedron_t(
    vertex_t const& v1,
    vertex_t const& v2,
    vertex_t const& v3,
    vertex_t const& v4)
    : v_{v1.vi(), v2.vi(), v3.vi(), v4.vi()}, edges_{}, faces_{}
{
}

std::array<index_type, 4u> const& tetrahedron_t::vertex_indices() const
{
    return v_;
}

std::array<index_type, 4u>& tetrahedron_t::vertex_indices()
{
    return v_;
}

index_type tetrahedron_t::v1() const
{
    return v_[0];
}

index_type tetrahedron_t::v2() const
{
    return v_[1];
}

index_type tetrahedron_t::v3() const
{
    return v_[2];
}

index_type tetrahedron_t::v4() const
{
    return v_[3];
}

index_type& tetrahedron_t::v1()
{
    return v_[0];
}

index_type& tetrahedron_t::v2()
{
    return v_[1];
}

index_type& tetrahedron_t::v3()
{
    return v_[2];
}

index_type& tetrahedron_t::v4()
{
    return v_[3];
}

std::array<edge_t, 6u> tetrahedron_t::edges_copy() const
{
    return std::array<edge_t, 6u>{
        edge_t{v1(), v2()},
        edge_t{v2(), v3()},
        edge_t{v3(), v1()},
        edge_t{v1(), v4()},
        edge_t{v2(), v4()},
        edge_t{v3(), v4()}};
}

std::array<triangle_t, 4u> tetrahedron_t::faces_copy() const
{
    return std::array<triangle_t, 4u>{
        triangle_t{v1(), v2(), v4()},
        triangle_t{v2(), v3(), v4()},
        triangle_t{v3(), v1(), v4()},
        triangle_t{v1(), v3(), v2()}};
}

std::array<index_type, 6u> const& tetrahedron_t::edge_indices() const
{
    return edges_;
}

std::array<index_type, 4u> const& tetrahedron_t::face_indices() const
{
    return faces_;
}

std::array<index_type, 6u>& tetrahedron_t::edge_indices()
{
    return edges_;
}

std::array<index_type, 4u>& tetrahedron_t::face_indices()
{
    return faces_;
}

std::uint8_t tetrahedron_t::id_of_vertex(index_type vi) const
{
    return static_cast<std::uint8_t>(
        std::distance(v_.begin(), std::find(v_.begin(), v_.end(), vi)));
}

std::uint8_t tetrahedron_t::id_of_edge(index_type ei) const
{
    return static_cast<std::uint8_t>(
        std::distance(edges_.begin(), std::find(edges_.begin(), edges_.end(), ei)));
}

std::uint8_t tetrahedron_t::id_of_face(index_type fi) const
{
    return static_cast<std::uint8_t>(
        std::distance(faces_.begin(), std::find(faces_.begin(), faces_.end(), fi)));
}

bool tetrahedron_t::operator==(tetrahedron_t const& other) const
{
    std::array<index_type, 4u> v      = v_;
    std::array<index_type, 4u> otherv = other.v_;

    std::sort(v.begin(), v.end());
    std::sort(otherv.begin(), otherv.end());

    return v == otherv;
}

bool tetrahedron_t::operator<(tetrahedron_t const& other) const
{
    std::array<index_type, 4u> v      = v_;
    std::array<index_type, 4u> otherv = other.v_;

    std::sort(v.begin(), v.end());
    std::sort(otherv.begin(), otherv.end());

    return v < otherv;
}

/**
 * Vertex set implementation
 */

vertex_t const& vertex_set_t::vertex(std::size_t vi) const
{
    return vertices_[vi];
}

vertex_t& vertex_set_t::vertex(std::size_t vi)
{
    return vertices_[vi];
}

bool vertex_set_t::contains_vertex(vertex_t const& v) const
{
    return v.vi() < vertex_count();
}

index_type vertex_set_t::add_vertex()
{
    index_type const vi = static_cast<index_type>(vertices_.size());
    vertices_.push_back(vertex_t{vi});

    return vi;
}

void vertex_set_t::add_vertex(vertex_t const& vertex)
{
    if (contains_vertex(vertex))
        return;

    std::size_t const previous_size = vertex_count();
    std::size_t const new_size      = static_cast<std::size_t>(vertex.vi()) + 1u;
    reserve_vertices(new_size);
    for (std::size_t i = previous_size; i < new_size; ++i)
    {
        add_vertex();
    }
}

void vertex_set_t::add_vertex(index_type const vi)
{
    add_vertex(vertex_t{vi});
}

std::size_t vertex_set_t::vertex_count() const
{
    return vertices_.size();
}

void vertex_set_t::reserve_vertices(std::size_t count)
{
    vertices_.reserve(count);
}

void vertex_set_t::clear()
{
    vertices_.clear();
}

std::vector<vertex_t> const& vertex_set_t::vertices() const
{
    return vertices_;
}

void vertex_set_t::remove_vertex_to_edge_incidency(index_type const vi, index_type ei)
{
    vertex_t& v = vertices_[vi];
    auto it = std::remove(v.incident_edge_indices().begin(), v.incident_edge_indices().end(), ei);
    v.incident_edge_indices().erase(it);
}

void vertex_set_t::remove_vertex_to_triangle_incidency(index_type const vi, index_type fi)
{
    vertex_t& v = vertices_[vi];
    auto it =
        std::remove(v.incident_triangle_indices().begin(), v.incident_triangle_indices().end(), fi);
    v.incident_triangle_indices().erase(it);
}

void vertex_set_t::remove_vertex_to_tetrahedron_incidency(index_type const vi, index_type ti)
{
    vertex_t& v = vertices_[vi];
    auto it     = std::remove(
        v.incident_tetrahedron_indices().begin(),
        v.incident_tetrahedron_indices().end(),
        ti);
    v.incident_tetrahedron_indices().erase(it);
}

bool vertex_set_t::operator==(vertex_set_t const& other) const
{
    return vertex_count() == other.vertex_count();
}

/**
 * Edge set implementation
 */

edge_t const& edge_set_t::edge(index_type ei) const
{
    return edges_[ei];
}

edge_t& edge_set_t::edge(index_type ei)
{
    return edges_[ei];
}

index_type edge_set_t::ei(edge_t const& edge) const
{
    auto it = edge_map_.find(edge);
    assert(it != edge_map_.end());
    return it->second;
}

index_type edge_set_t::add_edge(edge_t const& edge)
{
    auto it = edge_map_.find(edge);
    if (it != edge_map_.end())
    {
        index_type const ei = it->second;
        return ei;
    }

    edge_t created_edge{edge.v1(), edge.v2()};

    for (index_type const vi : created_edge.vertex_indices())
    {
        add_vertex(vi);
    }

    index_type ei{};
    if (!edge_garbage_collector_.empty())
    {
        ei         = *edge_garbage_collector_.begin();
        edges_[ei] = created_edge;
        edge_garbage_collector_.erase(edge_garbage_collector_.begin());
    }
    else
    {
        ei = static_cast<index_type>(edges_.size());
        edges_.push_back(created_edge);
    }

    edge_map_[created_edge] = ei;
    create_vertex_to_edge_incidency(ei);
    return ei;
}

edge_t edge_set_t::remove_edge(index_type ei)
{
    edge_t const& edge = edges_[ei];

    for (index_type const vi : edge.vertex_indices())
    {
        remove_vertex_to_edge_incidency(vi, ei);
    }

    auto it = edge_map_.find(edge);
    assert(it != edge_map_.end());
    edge_map_.erase(it);

    edge_garbage_collector_.insert(ei);

    return edge;
}

edge_t edge_set_t::remove_edge(edge_t const& edge)
{
    auto it = edge_map_.find(edge);
    assert(it != edge_map_.end());
    index_type const ei = it->second;
    remove_edge(ei);
    return edge;
}

std::size_t edge_set_t::edge_count() const
{
    return edges_.size() - edge_garbage_collector_.size();
}

bool edge_set_t::contains_edge(edge_t const& edge) const
{
    return edge_map_.find(edge) != edge_map_.end();
}

void edge_set_t::reserve_edges(std::size_t count)
{
    edges_.reserve(count);
}

void edge_set_t::clear()
{
    vertex_set_t::clear();
    edges_.clear();
    edge_map_.clear();
    edge_garbage_collector_.clear();
}

bool edge_set_t::is_safe_to_iterate_over_edges() const
{
    return edge_garbage_collector_.empty();
}

std::vector<edge_t> const& edge_set_t::edges() const
{
    return edges_;
}

std::vector<edge_t>& edge_set_t::edges()
{
    return edges_;
}

std::map<edge_t, index_type>::const_iterator edge_set_t::safe_edges_begin() const
{
    return edge_map_.begin();
}

std::map<edge_t, index_type>::const_iterator edge_set_t::safe_edges_end() const
{
    return edge_map_.end();
}

void edge_set_t::remove_edge_to_triangle_incidency(index_type const ei, index_type fi)
{
    edge_t& edge = edges_[ei];
    auto it      = std::remove(
        edge.incident_triangle_indices().begin(),
        edge.incident_triangle_indices().end(),
        fi);
    edge.incident_triangle_indices().erase(it);
}

void edge_set_t::remove_edge_to_tetrahedron_incidency(index_type const ei, index_type ti)
{
    edge_t& edge = edges_[ei];
    auto it      = std::remove(
        edge.incident_tetrahedron_indices().begin(),
        edge.incident_tetrahedron_indices().end(),
        ti);
    edge.incident_tetrahedron_indices().erase(it);
}

void edge_set_t::create_vertex_to_edge_incidency(index_type const ei)
{
    edge_t const& edge = edges_[ei];
    for (index_type const vi : edge.vertex_indices())
    {
        vertex_t& v = vertex(vi);
        v.incident_edge_indices().push_back(ei);
    }
}

bool edge_set_t::operator==(edge_set_t const& other) const
{
    bool const are_vertex_sets_equal = vertex_set_t::operator==(other);

    if (!are_vertex_sets_equal)
        return false;

    std::vector<edge_t> edge_set_1{};
    std::transform(
        edge_map_.begin(),
        edge_map_.end(),
        std::back_inserter(edge_set_1),
        [](std::pair<edge_t const, index_type> const& kv) { return kv.first; });

    std::vector<edge_t> edge_set_2{};
    std::transform(
        other.edge_map_.begin(),
        other.edge_map_.end(),
        std::back_inserter(edge_set_2),
        [](std::pair<edge_t const, index_type> const& kv) { return kv.first; });

    bool const are_edge_sets_equal = (edge_set_1 == edge_set_2);
    return are_edge_sets_equal;
}

/**
 * Triangle set implementation
 */

triangle_t const& triangle_set_t::triangle(index_type fi) const
{
    return triangles_[fi];
}

triangle_t& triangle_set_t::triangle(index_type fi)
{
    return triangles_[fi];
}

index_type triangle_set_t::fi(triangle_t const& triangle) const
{
    auto it = triangle_map_.find(triangle);
    assert(it != triangle_map_.end());
    return it->second;
}

void triangle_set_t::add_triangle(triangle_t const& triangle)
{
    auto it = triangle_map_.find(triangle);
    if (it != triangle_map_.end())
    {
        return;
    }

    triangle_t created_triangle{triangle.v1(), triangle.v2(), triangle.v3()};
    std::array<edge_t, 3u> const edge_copies = created_triangle.edges_copy();
    for (std::uint8_t e = 0u; e < edge_copies.size(); ++e)
    {
        edge_t const& edge_copy = edge_copies[e];

        auto edge_it = edge_map_.find(edge_copy);
        if (edge_it == edge_map_.end())
        {
            index_type const ei                   = add_edge(edge_copy);
            created_triangle.edge_indices().at(e) = ei;
        }
        else
        {
            index_type const ei                   = edge_it->second;
            created_triangle.edge_indices().at(e) = ei;
        }
    }
    index_type fi{};
    if (!triangle_garbage_collector_.empty())
    {
        fi             = *triangle_garbage_collector_.begin();
        triangles_[fi] = created_triangle;
        triangle_garbage_collector_.erase(triangle_garbage_collector_.begin());
    }
    else
    {
        fi = static_cast<index_type>(triangles_.size());
        triangles_.push_back(created_triangle);
    }

    triangle_map_[created_triangle] = fi;
    create_vertex_to_triangle_incidency(fi);
    create_edge_to_triangle_incidency(fi);
    return;
}

triangle_t triangle_set_t::remove_triangle(index_type fi)
{
    triangle_t const triangle = triangles_[fi];

    for (index_type const vi : triangle.vertex_indices())
    {
        remove_vertex_to_triangle_incidency(vi, fi);
    }
    for (index_type const ei : triangle.edge_indices())
    {
        remove_edge_to_triangle_incidency(ei, fi);
        if (edge(ei).incident_triangle_indices().empty())
        {
            remove_edge(ei);
        }
    }

    auto it = triangle_map_.find(triangle);
    assert(it != triangle_map_.end());
    triangle_map_.erase(it);
    triangle_garbage_collector_.insert(fi);

    return triangle;
}

triangle_t triangle_set_t::remove_triangle(triangle_t const& triangle)
{
    auto it = triangle_map_.find(triangle);
    assert(it != triangle_map_.end());
    index_type const fi = it->second;
    remove_triangle(fi);
    return triangle;
}

std::size_t triangle_set_t::triangle_count() const
{
    return triangles_.size() - triangle_garbage_collector_.size();
}

bool triangle_set_t::contains_triangle(triangle_t const& triangle) const
{
    return triangle_map_.find(triangle) != triangle_map_.end();
}

void triangle_set_t::reserve_triangles(std::size_t count)
{
    triangles_.reserve(count);
}

void triangle_set_t::clear()
{
    edge_set_t::clear();
    triangles_.clear();
    triangle_map_.clear();
    triangle_garbage_collector_.clear();
}

bool triangle_set_t::is_safe_to_iterate_over_triangles() const
{
    return triangle_garbage_collector_.empty();
}

std::vector<triangle_t> const& triangle_set_t::triangles() const
{
    return triangles_;
}

std::vector<triangle_t>& triangle_set_t::triangles()
{
    return triangles_;
}

std::map<triangle_t, index_type>::const_iterator triangle_set_t::safe_triangles_begin() const
{
    return triangle_map_.begin();
}

std::map<triangle_t, index_type>::const_iterator triangle_set_t::safe_triangles_end() const
{
    return triangle_map_.end();
}

void triangle_set_t::create_vertex_to_triangle_incidency(index_type const fi)
{
    triangle_t const& f = triangle(fi);
    for (index_type const vi : f.vertex_indices())
    {
        vertex_t& v = vertex(vi);
        v.incident_triangle_indices().push_back(fi);
    }
}

void triangle_set_t::create_edge_to_triangle_incidency(index_type const fi)
{
    triangle_t const& f = triangle(fi);
    for (index_type const ei : f.edge_indices())
    {
        edge_t& e = edge(ei);
        e.incident_triangle_indices().push_back(fi);
    }
}

void triangle_set_t::remove_triangle_to_tetrahedron_incidency(
    index_type const fi,
    index_type const ti)
{
    triangle_t& f = triangle(fi);
    auto it       = std::remove(
        f.incident_tetrahedron_indices().begin(),
        f.incident_tetrahedron_indices().end(),
        ti);
    f.incident_tetrahedron_indices().erase(it);
}

bool triangle_set_t::operator==(triangle_set_t const& other) const
{
    bool const are_edge_sets_equal = edge_set_t::operator==(other);

    if (!are_edge_sets_equal)
        return false;

    std::vector<triangle_t> triangle_set_1{};
    std::transform(
        triangle_map_.begin(),
        triangle_map_.end(),
        std::back_inserter(triangle_set_1),
        [](std::pair<triangle_t const, index_type> const& kv) { return kv.first; });

    std::vector<triangle_t> triangle_set_2{};
    std::transform(
        other.triangle_map_.begin(),
        other.triangle_map_.end(),
        std::back_inserter(triangle_set_2),
        [](std::pair<triangle_t const, index_type> const& kv) { return kv.first; });

    bool const are_triangle_sets_equal = (triangle_set_1 == triangle_set_2);
    return are_triangle_sets_equal;
}

/**
 * Tetrahedron set
 */

tetrahedron_t const& tetrahedron_set_t::tetrahedron(index_type ti) const
{
    return tetrahedra_[ti];
}

tetrahedron_t& tetrahedron_set_t::tetrahedron(index_type ti)
{
    return tetrahedra_[ti];
}

void tetrahedron_set_t::add_tetrahedron(tetrahedron_t const& tetrahedron)
{
    tetrahedron_t created_tetrahedron{
        tetrahedron.v1(),
        tetrahedron.v2(),
        tetrahedron.v3(),
        tetrahedron.v4()};

    std::array<triangle_t, 4u> const triangle_copies = created_tetrahedron.faces_copy();
    for (std::uint8_t f = 0u; f < triangle_copies.size(); ++f)
    {
        triangle_t const& triangle_copy = triangle_copies[f];

        auto triangle_it = triangle_map_.find(triangle_copy);
        if (triangle_it == triangle_map_.end())
        {
            add_triangle(triangle_copy);
            index_type const fi                      = triangle_map_[triangle_copy];
            created_tetrahedron.face_indices().at(f) = fi;
        }
        else
        {
            index_type const fi                      = triangle_it->second;
            created_tetrahedron.face_indices().at(f) = fi;
        }
    }
    std::array<edge_t, 6u> const edge_copies = created_tetrahedron.edges_copy();
    for (std::uint8_t e = 0u; e < edge_copies.size(); ++e)
    {
        auto it = edge_map_.find(edge_copies[e]);
        assert(it != edge_map_.end());
        index_type const ei                      = it->second;
        created_tetrahedron.edge_indices().at(e) = ei;
    }

    index_type ti{};
    if (!tetrahedron_garbage_collector_.empty())
    {
        ti              = *tetrahedron_garbage_collector_.begin();
        tetrahedra_[ti] = created_tetrahedron;
        tetrahedron_garbage_collector_.erase(tetrahedron_garbage_collector_.begin());
    }
    else
    {
        ti = static_cast<index_type>(tetrahedra_.size());
        tetrahedra_.push_back(created_tetrahedron);
    }

    create_vertex_to_tetrahedron_incidency(ti);
    create_edge_to_tetrahedron_incidency(ti);
    create_triangle_to_tetrahedron_incidency(ti);
    return;
}

tetrahedron_t tetrahedron_set_t::remove_tetrahedron(index_type ti)
{
    tetrahedron_t const tetrahedron = tetrahedra_[ti];

    for (index_type const vi : tetrahedron.vertex_indices())
    {
        remove_vertex_to_tetrahedron_incidency(vi, ti);
    }
    for (index_type const ei : tetrahedron.edge_indices())
    {
        remove_edge_to_tetrahedron_incidency(ei, ti);
    }
    for (index_type const fi : tetrahedron.face_indices())
    {
        remove_triangle_to_tetrahedron_incidency(fi, ti);
        if (triangle(fi).incident_tetrahedron_indices().empty())
        {
            remove_triangle(fi);
        }
    }

    tetrahedron_garbage_collector_.insert(ti);
    return tetrahedron;
}

std::size_t tetrahedron_set_t::tetrahedron_count() const
{
    return tetrahedra_.size() - tetrahedron_garbage_collector_.size();
}

void tetrahedron_set_t::reserve_tetrahedra(std::size_t count)
{
    tetrahedra_.reserve(count);
}

void tetrahedron_set_t::clear()
{
    triangle_set_t::clear();
    tetrahedra_.clear();
    tetrahedron_garbage_collector_.clear();
}

bool tetrahedron_set_t::is_safe_to_iterate_over_tetrahedra() const
{
    return tetrahedron_garbage_collector_.empty();
}

std::vector<tetrahedron_t> const& tetrahedron_set_t::tetrahedra() const
{
    return tetrahedra_;
}

std::vector<tetrahedron_t>& tetrahedron_set_t::tetrahedra()
{
    return tetrahedra_;
}

void tetrahedron_set_t::create_vertex_to_tetrahedron_incidency(index_type const ti)
{
    tetrahedron_t const& t = tetrahedron(ti);
    for (index_type const vi : t.vertex_indices())
    {
        vertex_t& v = vertex(vi);
        v.incident_tetrahedron_indices().push_back(ti);
    }
}

void tetrahedron_set_t::create_edge_to_tetrahedron_incidency(index_type const ti)
{
    tetrahedron_t const& t = tetrahedron(ti);
    for (index_type const ei : t.edge_indices())
    {
        edge_t& e = edge(ei);
        e.incident_tetrahedron_indices().push_back(ti);
    }
}

void tetrahedron_set_t::create_triangle_to_tetrahedron_incidency(index_type const ti)
{
    tetrahedron_t const& t = tetrahedron(ti);
    for (index_type const fi : t.face_indices())
    {
        triangle_t& f = triangle(fi);
        f.incident_tetrahedron_indices().push_back(ti);
    }
}

void tetrahedron_set_t::collect_garbage()
{
    /**
     * Partitions the edge vector, triangle vector and tetrahedron
     * vector such that primitives to remove are found at the end
     * of the vectors, and we can call erase() from the cutoff to
     * the end of the primitive vectors to delete those primitives.
     */

    std::size_t const previous_edge_count = edges().size();
    std::size_t const edge_cutoff         = edges().size() - edge_garbage_collector_.size();
    if (edge_cutoff > 0u)
    {
        auto pointer_to_valid_edge = edges().rbegin();
        auto const begin           = edges().begin();
        auto const index_of        = [this, begin](std::vector<edge_t>::reverse_iterator rit) {
            return static_cast<index_type>(std::distance(begin, rit.base())) - 1u;
        };
        auto const is_valid_edge = [this](index_type const ei) {
            return edge_garbage_collector_.find(ei) == edge_garbage_collector_.end();
        };
        // move the pointer to the first valid edge
        while (!is_valid_edge(index_of(pointer_to_valid_edge)))
        {
            ++pointer_to_valid_edge;
        }

        auto pointer_to_removed_edge = std::find_if(
            edge_garbage_collector_.rbegin(),
            edge_garbage_collector_.rend(),
            [pointer_to_valid_edge, index_of](index_type const ei) {
                index_type const eip = index_of(pointer_to_valid_edge);
                return ei < eip;
            });

        while (pointer_to_removed_edge != edge_garbage_collector_.rend())
        {
            index_type const ei  = *pointer_to_removed_edge;
            index_type const eip = index_of(pointer_to_valid_edge);

            swap_edges(ei, eip);

            ++pointer_to_valid_edge;
            ++pointer_to_removed_edge;
        }
    }

    std::size_t const previous_triangle_count = triangles().size();
    std::size_t const triangle_cutoff = triangles().size() - triangle_garbage_collector_.size();
    if (triangle_cutoff > 0u)
    {
        auto pointer_to_valid_triangle = triangles().rbegin();
        auto const begin               = triangles().begin();
        auto const index_of = [this, begin](std::vector<triangle_t>::reverse_iterator rit) {
            return static_cast<index_type>(std::distance(begin, rit.base())) - 1u;
        };
        auto const is_valid_triangle = [this](index_type const fi) {
            return triangle_garbage_collector_.find(fi) == triangle_garbage_collector_.end();
        };
        // move the pointer to the first valid edge
        while (!is_valid_triangle(index_of(pointer_to_valid_triangle)))
        {
            ++pointer_to_valid_triangle;
        }

        auto pointer_to_removed_triangle = std::find_if(
            triangle_garbage_collector_.rbegin(),
            triangle_garbage_collector_.rend(),
            [pointer_to_valid_triangle, index_of](index_type const fi) {
                index_type const fip = index_of(pointer_to_valid_triangle);
                return fi < fip;
            });

        while (pointer_to_removed_triangle != triangle_garbage_collector_.rend())
        {
            index_type const fi  = *pointer_to_removed_triangle;
            index_type const fip = index_of(pointer_to_valid_triangle);

            swap_triangles(fi, fip);

            ++pointer_to_valid_triangle;
            ++pointer_to_removed_triangle;
        }
    }

    std::size_t const previous_tetrahedron_count = tetrahedra().size();
    std::size_t const tetrahedron_cutoff =
        tetrahedra().size() - tetrahedron_garbage_collector_.size();
    if (tetrahedron_cutoff > 0u)
    {
        auto pointer_to_valid_tetrahedron = tetrahedra().rbegin();
        auto const begin                  = tetrahedra().begin();
        auto const index_of = [this, begin](std::vector<tetrahedron_t>::reverse_iterator rit) {
            return static_cast<index_type>(std::distance(begin, rit.base())) - 1u;
        };
        auto const is_valid_tetrahedron = [this](index_type const ti) {
            return tetrahedron_garbage_collector_.find(ti) == tetrahedron_garbage_collector_.end();
        };
        // move the pointer to the first valid edge
        while (!is_valid_tetrahedron(index_of(pointer_to_valid_tetrahedron)))
        {
            ++pointer_to_valid_tetrahedron;
        }

        auto pointer_to_removed_tetrahedron = std::find_if(
            tetrahedron_garbage_collector_.rbegin(),
            tetrahedron_garbage_collector_.rend(),
            [pointer_to_valid_tetrahedron, index_of](index_type const ti) {
                index_type const tip = index_of(pointer_to_valid_tetrahedron);
                return ti < tip;
            });

        while (pointer_to_removed_tetrahedron != tetrahedron_garbage_collector_.rend() &&
               pointer_to_valid_tetrahedron != tetrahedra().rend())
        {
            index_type const ti  = *pointer_to_removed_tetrahedron;
            index_type const tip = index_of(pointer_to_valid_tetrahedron);

            swap_tetrahedra(ti, tip);

            ++pointer_to_valid_tetrahedron;
            ++pointer_to_removed_tetrahedron;
        }
    }

    edges().erase(edges().begin() + edge_cutoff, edges().end());
    triangles().erase(triangles().begin() + triangle_cutoff, triangles().end());
    tetrahedra().erase(tetrahedra().begin() + tetrahedron_cutoff, tetrahedra().end());

    edge_garbage_collector_.clear();
    triangle_garbage_collector_.clear();
    tetrahedron_garbage_collector_.clear();
    return;
}

bool tetrahedron_set_t::is_boundary_tetrahedron(index_type const ti) const
{
    tetrahedron_t const& t = tetrahedron(ti);
    return is_boundary_tetrahedron(t);
}

bool tetrahedron_set_t::is_boundary_triangle(index_type const fi) const
{
    triangle_t const& f = triangle(fi);
    return is_boundary_triangle(f);
}

bool tetrahedron_set_t::is_boundary_edge(index_type const ei) const
{
    edge_t const& e = edge(ei);
    return is_boundary_edge(e);
}

bool tetrahedron_set_t::is_boundary_vertex(index_type const vi) const
{
    vertex_t const& v = vertex(vi);
    return is_boundary_vertex(v);
}

bool tetrahedron_set_t::is_boundary_tetrahedron(tetrahedron_t const& t) const
{
    auto const num_boundary_vertices = std::find_if(
        t.vertex_indices().begin(),
        t.vertex_indices().end(),
        [this](index_type const vi) { return is_boundary_vertex(vi); });
    return num_boundary_vertices != t.vertex_indices().end();
}

bool tetrahedron_set_t::is_boundary_triangle(triangle_t const& f) const
{
    return f.incident_tetrahedron_indices().size() == 1u;
}

bool tetrahedron_set_t::is_boundary_edge(edge_t const& e) const
{
    auto const num_boundary_triangles = std::count_if(
        e.incident_triangle_indices().begin(),
        e.incident_triangle_indices().end(),
        [this](index_type const fi) { return is_boundary_triangle(fi); });
    return num_boundary_triangles == 2u;
}

bool tetrahedron_set_t::is_boundary_vertex(vertex_t const& v) const
{
    auto it = std::find_if(
        v.incident_triangle_indices().begin(),
        v.incident_triangle_indices().end(),
        [this](index_type const fi) { return is_boundary_triangle(fi); });
    return it != v.incident_triangle_indices().end();
}

std::vector<tetrahedron_t> tetrahedron_set_t::boundary_tetrahedra() const
{
    std::vector<tetrahedron_t> boundary_tets{};
    for (tetrahedron_t const& t : tetrahedra_)
    {
        if (is_boundary_tetrahedron(t))
            boundary_tets.push_back(t);
    }
    return boundary_tets;
}

std::vector<index_type> tetrahedron_set_t::boundary_tetrahedron_indices() const
{
    std::vector<index_type> boundary_tet_indices{};
    for (std::size_t i = 0u; i < tetrahedra_.size(); ++i)
    {
        index_type const ti = static_cast<index_type>(i);
        if (is_boundary_tetrahedron(ti))
            boundary_tet_indices.push_back(ti);
    }
    return boundary_tet_indices;
}

std::vector<triangle_t> tetrahedron_set_t::boundary_triangles() const
{
    std::vector<triangle_t> boundary_tris{};
    for (triangle_t const& f : triangles())
    {
        if (is_boundary_triangle(f))
            boundary_tris.push_back(f);
    }
    return boundary_tris;
}

std::vector<index_type> tetrahedron_set_t::boundary_triangle_indices() const
{
    std::vector<index_type> boundary_triangle_indices{};
    for (std::size_t i = 0u; i < triangles().size(); ++i)
    {
        index_type const fi = static_cast<index_type>(i);
        if (is_boundary_triangle(fi))
            boundary_triangle_indices.push_back(fi);
    }
    return boundary_triangle_indices;
}

std::vector<edge_t> tetrahedron_set_t::boundary_edges() const
{
    std::vector<edge_t> boundary_edge_list{};
    for (edge_t const& e : edges())
    {
        if (is_boundary_edge(e))
            boundary_edge_list.push_back(e);
    }
    return boundary_edge_list;
}

std::vector<index_type> tetrahedron_set_t::boundary_edge_indices() const
{
    std::vector<index_type> boundary_edge_indices{};
    for (std::size_t i = 0u; i < edges().size(); ++i)
    {
        index_type const ei = static_cast<index_type>(i);
        if (is_boundary_edge(ei))
            boundary_edge_indices.push_back(ei);
    }
    return boundary_edge_indices;
}

std::vector<vertex_t> tetrahedron_set_t::boundary_vertices() const
{
    std::vector<vertex_t> boundary_vertex_list{};
    for (vertex_t const& v : vertices())
    {
        if (is_boundary_vertex(v))
            boundary_vertex_list.push_back(v);
    }
    return boundary_vertex_list;
}

std::vector<index_type> tetrahedron_set_t::boundary_vertex_indices() const
{
    std::vector<index_type> boundary_vertex_indices{};
    for (std::size_t i = 0u; i < vertices().size(); ++i)
    {
        index_type const vi = static_cast<index_type>(i);
        if (is_boundary_vertex(vi))
            boundary_vertex_indices.push_back(vi);
    }
    return boundary_vertex_indices;
}

std::vector<tetrahedron_t> tetrahedron_set_t::interior_tetrahedra() const
{
    std::vector<tetrahedron_t> interior_tets{};
    for (tetrahedron_t const& t : tetrahedra_)
    {
        if (!is_boundary_tetrahedron(t))
            interior_tets.push_back(t);
    }
    return interior_tets;
}

std::vector<index_type> tetrahedron_set_t::interior_tetrahedron_indices() const
{
    std::vector<index_type> interior_tet_indices{};
    for (std::size_t i = 0u; i < tetrahedra_.size(); ++i)
    {
        index_type const ti = static_cast<index_type>(i);
        if (!is_boundary_tetrahedron(ti))
            interior_tet_indices.push_back(ti);
    }
    return interior_tet_indices;
}

std::vector<triangle_t> tetrahedron_set_t::interior_triangles() const
{
    std::vector<triangle_t> interior_tris{};
    for (triangle_t const& f : triangles())
    {
        if (!is_boundary_triangle(f))
            interior_tris.push_back(f);
    }
    return interior_tris;
}

std::vector<index_type> tetrahedron_set_t::interior_triangle_indices() const
{
    std::vector<index_type> interior_triangle_indices{};
    for (std::size_t i = 0u; i < triangles().size(); ++i)
    {
        index_type const fi = static_cast<index_type>(i);
        if (!is_boundary_triangle(fi))
            interior_triangle_indices.push_back(fi);
    }
    return interior_triangle_indices;
}

std::vector<edge_t> tetrahedron_set_t::interior_edges() const
{
    std::vector<edge_t> interior_edge_list{};
    for (edge_t const& e : edges())
    {
        if (!is_boundary_edge(e))
            interior_edge_list.push_back(e);
    }
    return interior_edge_list;
}

std::vector<index_type> tetrahedron_set_t::interior_edge_indices() const
{
    std::vector<index_type> interior_edge_indices{};
    for (std::size_t i = 0u; i < edges().size(); ++i)
    {
        index_type const ei = static_cast<index_type>(i);
        if (!is_boundary_edge(ei))
            interior_edge_indices.push_back(ei);
    }
    return interior_edge_indices;
}

std::vector<vertex_t> tetrahedron_set_t::interior_vertices() const
{
    std::vector<vertex_t> interior_vertex_list{};
    for (vertex_t const& v : vertices())
    {
        if (!is_boundary_vertex(v))
            interior_vertex_list.push_back(v);
    }
    return interior_vertex_list;
}

std::vector<index_type> tetrahedron_set_t::interior_vertex_indices() const
{
    std::vector<index_type> interior_vertex_indices{};
    for (std::size_t i = 0u; i < vertices().size(); ++i)
    {
        index_type const vi = static_cast<index_type>(i);
        if (!is_boundary_vertex(vi))
            interior_vertex_indices.push_back(vi);
    }
    return interior_vertex_indices;
}

std::vector<triangle_t> tetrahedron_set_t::oriented_boundary_triangles() const
{
    std::vector<triangle_t> oriented_boundary_tris{};
    for (std::size_t i = 0u; i < triangles().size(); ++i)
    {
        index_type const fi = static_cast<index_type>(i);
        if (is_boundary_triangle(fi))
        {
            triangle_t const& f = triangle(fi);
            // boundary triangles only have one incident tet
            index_type const ti         = f.incident_tetrahedron_indices().front();
            tetrahedron_t const& t      = tetrahedron(ti);
            std::uint8_t const f_idx    = t.id_of_face(fi);
            triangle_t const f_oriented = t.faces_copy()[f_idx];
            oriented_boundary_tris.push_back(f_oriented);
        }
    }
    return oriented_boundary_tris;
}

bool tetrahedron_set_t::operator==(tetrahedron_set_t const& other) const
{
    bool const are_triangle_sets_equal = triangle_set_t::operator==(other);

    if (!are_triangle_sets_equal)
        return false;

    std::set<tetrahedron_t> tet_set_1(tetrahedra_.begin(), tetrahedra_.end());
    std::set<tetrahedron_t> tet_set_2(other.tetrahedra_.begin(), other.tetrahedra_.end());

    bool const are_tetrahedron_sets_equal = (tet_set_1 == tet_set_2);
    return are_tetrahedron_sets_equal;
}

void tetrahedron_set_t::swap_edges(index_type ei, index_type eip)
{
    std::swap(edge(ei), edge(eip));
    edge_map_[edge(ei)] = ei;

    edge_t const& swapped_edge = edge(ei);

    for (index_type const vi : swapped_edge.vertex_indices())
    {
        vertex_t& v = vertex(vi);
        std::replace(v.incident_edge_indices().begin(), v.incident_edge_indices().end(), eip, ei);
    }
    for (index_type const fi : swapped_edge.incident_triangle_indices())
    {
        triangle_t& f = triangle(fi);
        std::replace(f.edge_indices().begin(), f.edge_indices().end(), eip, ei);
    }
    for (index_type const ti : swapped_edge.incident_tetrahedron_indices())
    {
        tetrahedron_t& t = tetrahedron(ti);
        std::replace(t.edge_indices().begin(), t.edge_indices().end(), eip, ei);
    }
}

void tetrahedron_set_t::swap_triangles(index_type fi, index_type fip)
{
    std::swap(triangle(fi), triangle(fip));
    triangle_map_[triangle(fi)] = fi;

    triangle_t const& swapped_triangle = triangle(fi);

    for (index_type const vi : swapped_triangle.vertex_indices())
    {
        vertex_t& v = vertex(vi);
        std::replace(
            v.incident_triangle_indices().begin(),
            v.incident_triangle_indices().end(),
            fip,
            fi);
    }
    for (index_type const ei : swapped_triangle.edge_indices())
    {
        edge_t& e = edge(ei);
        std::replace(
            e.incident_triangle_indices().begin(),
            e.incident_triangle_indices().end(),
            fip,
            fi);
    }
    for (index_type const ti : swapped_triangle.incident_tetrahedron_indices())
    {
        tetrahedron_t& t = tetrahedron(ti);
        std::replace(t.face_indices().begin(), t.face_indices().end(), fip, fi);
    }
}

void tetrahedron_set_t::swap_tetrahedra(index_type ti, index_type tip)
{
    std::swap(tetrahedron(ti), tetrahedron(tip));

    tetrahedron_t const& swapped_tetrahedron = tetrahedron(ti);

    for (index_type const vi : swapped_tetrahedron.vertex_indices())
    {
        vertex_t& v = vertex(vi);
        std::replace(
            v.incident_tetrahedron_indices().begin(),
            v.incident_tetrahedron_indices().end(),
            tip,
            ti);
    }
    for (index_type const ei : swapped_tetrahedron.edge_indices())
    {
        edge_t& e = edge(ei);
        std::replace(
            e.incident_tetrahedron_indices().begin(),
            e.incident_tetrahedron_indices().end(),
            tip,
            ti);
    }
    for (index_type const fi : swapped_tetrahedron.face_indices())
    {
        triangle_t& f = triangle(fi);
        std::replace(
            f.incident_tetrahedron_indices().begin(),
            f.incident_tetrahedron_indices().end(),
            tip,
            ti);
    }
}

} // namespace topology
} // namespace sbs