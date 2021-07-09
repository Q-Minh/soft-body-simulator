#include "physics/mesh.h"

#include "common/geometry.h"
#include "common/primitive.h"

#include <map>
#include <numeric>

namespace sbs {
namespace physics {

vertex_t::vertex_t()
{
    p_.setZero();
    v_.setZero();
    f_.setZero();
    m_ = 0.;
    n_.setZero();
    c_.setZero();
}

vertex_t::position_type const& vertex_t::position() const
{
    return p_;
}

vertex_t::position_type& vertex_t::position()
{
    return p_;
}

vertex_t::velocity_type const& vertex_t::velocity() const
{
    return v_;
}

vertex_t::velocity_type& vertex_t::velocity()
{
    return v_;
}

vertex_t::force_type const& vertex_t::force() const
{
    return f_;
}

vertex_t::force_type& vertex_t::force()
{
    return f_;
}

scalar_type const& vertex_t::mass() const
{
    return m_;
}

scalar_type& vertex_t::mass()
{
    return m_;
}

vertex_t::normal_type const& vertex_t::normal() const
{
    return n_;
}

vertex_t::normal_type& vertex_t::normal()
{
    return n_;
}

vertex_t::color_type const& vertex_t::color() const
{
    return c_;
}

vertex_t::color_type& vertex_t::color()
{
    return c_;
}

tetrahedron_t::tetrahedron_t(index_type v1, index_type v2, index_type v3, index_type v4)
    : v_{v1, v2, v3, v4}, edges_{}, faces_{}, rho_{0.}
{
}

tetrahedron_t::tetrahedron_t(tetrahedron_t const& other)
    : v_{other.v_}, edges_{}, faces_{}, rho_{other.rho_}
{
    if (static_cast<bool>(other.edges_))
    {
        edges_ = std::make_unique<std::array<index_type, 6u>>(*other.edges_);
    }
    if (static_cast<bool>(other.faces_))
    {
        faces_ = std::make_unique<std::array<index_type, 4u>>(*other.faces_);
    }
}

index_type const& tetrahedron_t::v1() const
{
    return v_[0];
}

index_type const& tetrahedron_t::v2() const
{
    return v_[1];
}

index_type const& tetrahedron_t::v3() const
{
    return v_[2];
}

index_type const& tetrahedron_t::v4() const
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

scalar_type const& tetrahedron_t::mass_density() const
{
    return rho_;
}

scalar_type& tetrahedron_t::mass_density()
{
    return rho_;
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

std::unique_ptr<std::array<index_type, 6u>> const& tetrahedron_t::edge_indices() const
{
    return edges_;
}

std::unique_ptr<std::array<index_type, 6u>>& tetrahedron_t::edge_indices()
{
    return edges_;
}

std::unique_ptr<std::array<index_type, 4u>> const& tetrahedron_t::face_indices() const
{
    return faces_;
}

std::unique_ptr<std::array<index_type, 4u>>& tetrahedron_t::face_indices()
{
    return faces_;
}

edge_t::edge_t(index_type v1, index_type v2) : v_{v1, v2}, adjacent_tets_{}, adjacent_triangles_{}
{
}

index_type const& edge_t::v1() const
{
    return v_[0];
}

index_type const& edge_t::v2() const
{
    return v_[1];
}

index_type& edge_t::v1()
{
    return v_[0];
}

index_type& edge_t::v2()
{
    return v_[1];
}

std::vector<index_type> const& edge_t::adjacent_tetrahedron_indices() const
{
    return adjacent_tets_;
}

std::vector<index_type>& edge_t::adjacent_tetrahedron_indices()
{
    return adjacent_tets_;
}

std::vector<index_type> const& edge_t::adjacent_triangle_indices() const
{
    return adjacent_triangles_;
}

std::vector<index_type>& edge_t::adjacent_triangle_indices()
{
    return adjacent_triangles_;
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

bool edge_t::is_reverse_of(edge_t const& other) const
{
    return v1() == other.v2() && v2() == other.v1();
}

triangle_t::triangle_t(index_type v1, index_type v2, index_type v3)
    : v_{v1, v2, v3}, adjacent_tets_{}, edges_{}
{
}

triangle_t::triangle_t(triangle_t const& other)
    : v_{other.v_}, adjacent_tets_{other.adjacent_tets_}, edges_{}
{
    if (static_cast<bool>(other.edges_))
    {
        edges_ = std::make_unique<std::array<index_type, 3u>>(*other.edges_);
    }
}

index_type const& triangle_t::v1() const
{
    return v_[0];
}

index_type const& triangle_t::v2() const
{
    return v_[1];
}

index_type const& triangle_t::v3() const
{
    return v_[2];
}

index_type& triangle_t::v1()
{
    return v_[0];
}

index_type& triangle_t::v2()
{
    return v_[1];
}

index_type& triangle_t::v3()
{
    return v_[2];
}

std::array<edge_t, 3u> triangle_t::edges_copy() const
{
    return std::array<edge_t, 3u>{edge_t{v1(), v2()}, edge_t{v2(), v3()}, edge_t{v3(), v1()}};
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

bool triangle_t::is_boundary_triangle() const
{
    return adjacent_tets_.size() == 1u;
}

bool triangle_t::is_interior_triangle() const
{
    return adjacent_tets_.size() == 2u;
}

std::vector<index_type> const& triangle_t::adjacent_tetrahedron_indices() const
{
    return adjacent_tets_;
}

std::vector<index_type>& triangle_t::adjacent_tetrahedron_indices()
{
    return adjacent_tets_;
}

std::unique_ptr<std::array<index_type, 3u>> const& triangle_t::adjacent_edge_indices() const
{
    return edges_;
}

std::unique_ptr<std::array<index_type, 3u>>& triangle_t::adjacent_edge_indices()
{
    return edges_;
}

void build_topology_parameters_t::include_vertex_to_edge_adjacency()
{
    vertex_to_edge = true;
}
void build_topology_parameters_t::include_vertex_to_triangle_adjacency()
{
    vertex_to_triangle = true;
}
void build_topology_parameters_t::include_vertex_to_tetrahedra_adjacency()
{
    vertex_to_tetrahedra = true;
}
void build_topology_parameters_t::include_triangle_to_edge_adjacency()
{
    triangle_to_edge = true;
}
void build_topology_parameters_t::include_triangle_to_tetrahedra_adjacency()
{
    triangle_to_tetrahedra = true;
}
void build_topology_parameters_t::include_tetrahedron_to_edge_adjacency()
{
    tetrahedron_to_edge = true;
}
void build_topology_parameters_t::include_tetrahedron_to_triangle_adjacency()
{
    tetrahedron_to_triangle = true;
}
void build_topology_parameters_t::include_edge_to_triangle_adjacency()
{
    edge_to_triangle = true;
}
void build_topology_parameters_t::include_edge_to_tetrahedra_adjacency()
{
    edge_to_tetrahedra = true;
}

void build_topology_parameters_t::exclude_vertex_to_edge_adjacency()
{
    vertex_to_edge = false;
}
void build_topology_parameters_t::exclude_vertex_to_triangle_adjacency()
{
    vertex_to_triangle = false;
}
void build_topology_parameters_t::exclude_vertex_to_tetrahedra_adjacency()
{
    vertex_to_tetrahedra = false;
}
void build_topology_parameters_t::exclude_triangle_to_edge_adjacency()
{
    triangle_to_edge = false;
}
void build_topology_parameters_t::exclude_triangle_to_tetrahedra_adjacency()
{
    triangle_to_tetrahedra = false;
}
void build_topology_parameters_t::exclude_tetrahedron_to_edge_adjacency()
{
    tetrahedron_to_edge = false;
}
void build_topology_parameters_t::exclude_tetrahedron_to_triangle_adjacency()
{
    tetrahedron_to_triangle = false;
}
void build_topology_parameters_t::exclude_edge_to_triangle_adjacency()
{
    edge_to_triangle = false;
}
void build_topology_parameters_t::exclude_edge_to_tetrahedra_adjacency()
{
    edge_to_tetrahedra = false;
}

topological_simulated_tetrahedral_mesh_t::topological_simulated_tetrahedral_mesh_t(
    common::geometry_t const& geometry,
    build_topology_parameters_t const& params)
    : vertices_(), edges_(), triangles_(), tetrahedra_(), topology_params_(params)
{
    bool const is_tet_mesh =
        geometry.geometry_type == common::geometry_t::geometry_type_t::tetrahedron;
    assert(is_tet_mesh && "is tet mesh failed");
    if (!is_tet_mesh)
        return;

    assert(geometry.positions.size() % 3 == 0 && "must have (x,y,z) per position");
    assert(geometry.indices.size() % 4 == 0 && "must have (v1,v2,v3,v4) per tetrahedron");

    auto const vertex_count = geometry.positions.size() / 3u;
    vertices_.reserve(vertex_count);

    auto const tetrahedra_count = geometry.indices.size() / 4u;
    tetrahedra_.reserve(tetrahedra_count);

    bool const geometry_has_normals = geometry.has_normals();
    bool const geometry_has_colors  = geometry.has_colors();

    if (geometry_has_normals)
    {
        assert(
            geometry.normals.size() == geometry.positions.size() &&
            "must have as many normals as positions");
    }
    if (geometry_has_colors)
    {
        assert(
            geometry.colors.size() == geometry.positions.size() &&
            "must have as many colors as positions");
    }

    for (std::size_t i = 0u; i < geometry.positions.size(); i += 3u)
    {
        vertex_t vertex{};

        auto const x = geometry.positions[i];
        auto const y = geometry.positions[i + 1u];
        auto const z = geometry.positions[i + 2u];

        vertex.position() = vertex_t::position_type{x, y, z};

        if (geometry_has_normals)
        {
            auto const nx   = geometry.normals[i];
            auto const ny   = geometry.normals[i + 1u];
            auto const nz   = geometry.normals[i + 2u];
            vertex.normal() = vertex_t::normal_type{nx, ny, nz};
        }
        if (geometry_has_colors)
        {
            auto const r   = static_cast<float>(geometry.colors[i]) / 255.f;
            auto const g   = static_cast<float>(geometry.colors[i + 1u]) / 255.f;
            auto const b   = static_cast<float>(geometry.colors[i + 2u]) / 255.f;
            vertex.color() = vertex_t::color_type{r, g, b};
        }

        vertices_.push_back(vertex);
    }

    for (std::size_t i = 0u; i < geometry.indices.size(); i += 4u)
    {
        tetrahedron_t tetrahedron{};

        auto const v1 = geometry.indices[i];
        auto const v2 = geometry.indices[i + 1u];
        auto const v3 = geometry.indices[i + 2u];
        auto const v4 = geometry.indices[i + 3u];

        tetrahedron.v1() = v1;
        tetrahedron.v2() = v2;
        tetrahedron.v3() = v3;
        tetrahedron.v4() = v4;

        tetrahedra_.push_back(tetrahedron);
    }

    build_topology(topology_params_);
}

std::vector<vertex_t> const& topological_simulated_tetrahedral_mesh_t::vertices() const
{
    return vertices_;
}

std::vector<vertex_t>& topological_simulated_tetrahedral_mesh_t::vertices()
{
    return vertices_;
}

std::vector<edge_t> const& topological_simulated_tetrahedral_mesh_t::edges() const
{
    return edges_;
}

std::vector<edge_t>& topological_simulated_tetrahedral_mesh_t::edges()
{
    return edges_;
}

std::vector<triangle_t> const& topological_simulated_tetrahedral_mesh_t::triangles() const
{
    return triangles_;
}

std::vector<triangle_t>& topological_simulated_tetrahedral_mesh_t::triangles()
{
    return triangles_;
}

std::vector<tetrahedron_t> const& topological_simulated_tetrahedral_mesh_t::tetrahedra() const
{
    return tetrahedra_;
}

std::vector<tetrahedron_t>& topological_simulated_tetrahedral_mesh_t::tetrahedra()
{
    return tetrahedra_;
}

build_topology_parameters_t const&
topological_simulated_tetrahedral_mesh_t::topology_parameters() const
{
    return topology_params_;
}

void topological_simulated_tetrahedral_mesh_t::build_topology(
    build_topology_parameters_t const& params)
{
    edges_.clear();
    triangles_.clear();

    // pre-allocate memory to speed up construction. might overuse memory
    edges_.reserve(tetrahedra_.size() * 2u);
    triangles_.reserve(tetrahedra_.size() * 4u);

    std::map<edge_t, index_type /* index of edge in the edges_ vector */> edge_map{};
    std::map<triangle_t, index_type /* index of triangle in the triangles_ vector */>
        triangle_map{};

    // compute tetrahedron<->edge, tetrahedron<->triangle adjacencies
    for (std::size_t i = 0u; i < tetrahedra_.size(); ++i)
    {
        tetrahedron_t& tetrahedron = tetrahedra_[i];
        index_type const t         = static_cast<index_type>(i);

        tetrahedron.edge_indices().reset();
        tetrahedron.face_indices().reset();

        std::array<edge_t, 6u> const edge_copies     = tetrahedron.edges_copy();
        std::array<triangle_t, 4u> const face_copies = tetrahedron.faces_copy();

        if (params.tetrahedron_to_edge)
        {
            tetrahedron.edge_indices() = std::make_unique<std::array<index_type, 6u>>();
        }
        if (params.tetrahedron_to_triangle)
        {
            tetrahedron.face_indices() = std::make_unique<std::array<index_type, 4u>>();
        }

        for (std::size_t j = 0u; j < edge_copies.size(); ++j)
        {
            std::map<edge_t, index_type>::iterator edge_it = edge_map.find(edge_copies[j]);
            if (edge_it == edge_map.end())
            {
                std::size_t const e = edges_.size();
                edges_.push_back(edge_copies[j]);
                auto const [it, was_inserted] =
                    edge_map.insert(std::make_pair(edge_copies[j], static_cast<index_type>(e)));

                assert(was_inserted);

                edge_it = it;
            }

            index_type const e = edge_it->second;
            edge_t& edge       = edges_[e];

            if (params.edge_to_tetrahedra)
            {
                edge.adjacent_tetrahedron_indices().push_back(t);
            }
            if (params.tetrahedron_to_edge)
            {
                tetrahedron.edge_indices()->at(j) = e;
            }
        }

        for (std::size_t j = 0u; j < face_copies.size(); ++j)
        {
            std::map<triangle_t, index_type>::iterator face_it = triangle_map.find(face_copies[j]);
            if (face_it == triangle_map.end())
            {
                std::size_t const f = triangles_.size();
                triangles_.push_back(face_copies[j]);
                auto const [it, was_inserted] =
                    triangle_map.insert(std::make_pair(face_copies[j], static_cast<index_type>(f)));

                assert(was_inserted);
                face_it = it;
            }

            index_type const f   = face_it->second;
            triangle_t& triangle = triangles_[f];

            if (params.triangle_to_tetrahedra)
            {
                triangle.adjacent_tetrahedron_indices().push_back(t);
            }
            if (params.tetrahedron_to_triangle)
            {
                tetrahedron.face_indices()->at(j) = f;
            }
        }
    }

    // compute triangle<->edge adjacencies
    for (std::size_t i = 0u; i < triangles_.size(); ++i)
    {
        triangle_t& triangle = triangles_[i];
        index_type const f   = static_cast<index_type>(i);

        if (params.triangle_to_edge)
        {
            triangle.adjacent_edge_indices() = std::make_unique<std::array<index_type, 3u>>();

            std::array<edge_t, 3u> const edge_copies = triangle.edges_copy();

            for (std::size_t j = 0u; j < edge_copies.size(); ++j)
            {
                index_type const e                      = edge_map[edge_copies[j]];
                triangle.adjacent_edge_indices()->at(j) = e;
            }
        }
        if (params.edge_to_triangle)
        {
            std::array<edge_t, 3u> const edge_copies = triangle.edges_copy();

            for (std::size_t j = 0u; j < edge_copies.size(); ++j)
            {
                index_type const e = edge_map[edge_copies[j]];
                edges_[e].adjacent_triangle_indices().push_back(f);
            }
        }
    }
}

std::vector<triangle_t const*> topological_simulated_tetrahedral_mesh_t::boundary_triangles() const
{
    std::vector<triangle_t const*> triangles(triangles_.size());
    std::transform(
        triangles_.begin(),
        triangles_.end(),
        triangles.begin(),
        [](triangle_t const& t) { return &t; });

    std::size_t const num_boundary_triangles = boundary_triangle_count();
    std::vector<triangle_t const*> triangles_on_boundary(num_boundary_triangles);
    std::copy_if(
        triangles.begin(),
        triangles.end(),
        triangles_on_boundary.begin(),
        [](triangle_t const* t) { return t->is_boundary_triangle(); });

    return triangles_on_boundary;
}

std::size_t topological_simulated_tetrahedral_mesh_t::boundary_triangle_count() const
{
    return std::count_if(triangles().begin(), triangles().end(), [](triangle_t const& triangle) {
        return triangle.is_boundary_triangle();
    });
}

std::vector<tetrahedron_t const*>
topological_simulated_tetrahedral_mesh_t::edge_adjacent_tetrahedra(index_type e) const
{
    edge_t const& edge                                  = edges_[e];
    std::vector<index_type> const& adjacent_tet_indices = edge.adjacent_tetrahedron_indices();
    std::vector<tetrahedron_t const*> adjacent_tets(adjacent_tet_indices.size());
    std::transform(
        adjacent_tet_indices.begin(),
        adjacent_tet_indices.end(),
        adjacent_tets.begin(),
        [this](index_type const t) { return &tetrahedra_[t]; });
    return adjacent_tets;
}

std::vector<tetrahedron_t const*>
topological_simulated_tetrahedral_mesh_t::triangle_adjacent_tetrahedra(index_type f) const
{
    triangle_t const& triangle                          = triangles_[f];
    std::vector<index_type> const& adjacent_tet_indices = triangle.adjacent_tetrahedron_indices();
    std::vector<tetrahedron_t const*> adjacent_tets(adjacent_tet_indices.size());
    std::transform(
        adjacent_tet_indices.begin(),
        adjacent_tet_indices.end(),
        adjacent_tets.begin(),
        [this](index_type const t) { return &tetrahedra_[t]; });
    return adjacent_tets;
}

std::array<edge_t const*, 6u>
topological_simulated_tetrahedral_mesh_t::edges_of_tetrahedron(index_type t) const
{
    tetrahedron_t const& tetrahedron = tetrahedra_[t];
    assert(
        static_cast<bool>(tetrahedron.edge_indices()) &&
        "Must have tetrahedron to edge connectivity info");

    std::array<index_type, 6u> const& tet_edge_indices = *tetrahedron.edge_indices();
    std::array<edge_t const*, 6u> const tet_edges{
        &edges_[tet_edge_indices[0]],
        &edges_[tet_edge_indices[1]],
        &edges_[tet_edge_indices[2]],
        &edges_[tet_edge_indices[3]],
        &edges_[tet_edge_indices[4]],
        &edges_[tet_edge_indices[5]],
    };

    return tet_edges;
}

std::array<triangle_t const*, 4u>
topological_simulated_tetrahedral_mesh_t::faces_of_tetrahedron(index_type t) const
{
    tetrahedron_t const& tetrahedron = tetrahedra_[t];
    assert(
        static_cast<bool>(tetrahedron.face_indices()) &&
        "Must have tetrahedron to face connectivity info");

    std::array<index_type, 4u> const& tet_face_indices = *tetrahedron.face_indices();
    std::array<triangle_t const*, 4u> const tet_faces{
        &triangles_[tet_face_indices[0]],
        &triangles_[tet_face_indices[1]],
        &triangles_[tet_face_indices[2]],
        &triangles_[tet_face_indices[3]]};

    return tet_faces;
}

std::array<vertex_t const*, 4u>
topological_simulated_tetrahedral_mesh_t::vertices_of_tetrahedron(index_type t) const
{
    tetrahedron_t const& tetrahedron = tetrahedra_[t];
    return std::array<vertex_t const*, 4u>{
        &vertices_[tetrahedron.v1()],
        &vertices_[tetrahedron.v2()],
        &vertices_[tetrahedron.v3()],
        &vertices_[tetrahedron.v4()],
    };
}

std::vector<tetrahedron_t*>
topological_simulated_tetrahedral_mesh_t::edge_adjacent_tetrahedra(index_type e)
{
    edge_t const& edge                                  = edges_[e];
    std::vector<index_type> const& adjacent_tet_indices = edge.adjacent_tetrahedron_indices();
    std::vector<tetrahedron_t*> adjacent_tets(adjacent_tet_indices.size());
    std::transform(
        adjacent_tet_indices.begin(),
        adjacent_tet_indices.end(),
        adjacent_tets.begin(),
        [this](index_type const t) { return &tetrahedra_[t]; });
    return adjacent_tets;
}

std::vector<tetrahedron_t*>
topological_simulated_tetrahedral_mesh_t::triangle_adjacent_tetrahedra(index_type f)
{
    triangle_t const& triangle                          = triangles_[f];
    std::vector<index_type> const& adjacent_tet_indices = triangle.adjacent_tetrahedron_indices();
    std::vector<tetrahedron_t*> adjacent_tets(adjacent_tet_indices.size());
    std::transform(
        adjacent_tet_indices.begin(),
        adjacent_tet_indices.end(),
        adjacent_tets.begin(),
        [this](index_type const t) { return &tetrahedra_[t]; });
    return adjacent_tets;
}

std::array<edge_t*, 6u> topological_simulated_tetrahedral_mesh_t::edges_of_tetrahedron(index_type t)
{
    tetrahedron_t const& tetrahedron = tetrahedra_[t];
    assert(
        static_cast<bool>(tetrahedron.edge_indices()) &&
        "Must have tetrahedron to edge connectivity info");

    std::array<index_type, 6u> const& tet_edge_indices = *tetrahedron.edge_indices();
    std::array<edge_t*, 6u> const tet_edges{
        &edges_[tet_edge_indices[0]],
        &edges_[tet_edge_indices[1]],
        &edges_[tet_edge_indices[2]],
        &edges_[tet_edge_indices[3]],
        &edges_[tet_edge_indices[4]],
        &edges_[tet_edge_indices[5]],
    };

    return tet_edges;
}

std::array<triangle_t*, 4u>
topological_simulated_tetrahedral_mesh_t::faces_of_tetrahedron(index_type t)
{
    tetrahedron_t const& tetrahedron = tetrahedra_[t];
    assert(
        static_cast<bool>(tetrahedron.face_indices()) &&
        "Must have tetrahedron to face connectivity info");

    std::array<index_type, 4u> const& tet_face_indices = *tetrahedron.face_indices();
    std::array<triangle_t*, 4u> const tet_faces{
        &triangles_[tet_face_indices[0]],
        &triangles_[tet_face_indices[1]],
        &triangles_[tet_face_indices[2]],
        &triangles_[tet_face_indices[3]]};

    return tet_faces;
}

std::array<vertex_t*, 4u>
topological_simulated_tetrahedral_mesh_t::vertices_of_tetrahedron(index_type t)
{
    tetrahedron_t const& tetrahedron = tetrahedra_[t];
    return std::array<vertex_t*, 4u>{
        &vertices_[tetrahedron.v1()],
        &vertices_[tetrahedron.v2()],
        &vertices_[tetrahedron.v3()],
        &vertices_[tetrahedron.v4()],
    };
}

bool topological_simulated_tetrahedral_mesh_t::is_boundary_edge(index_type e) const
{
    // a boundary edge should have two adjacent boundary triangles
    edge_t const& edge = edges_[e];
    return std::count_if(
               edge.adjacent_triangle_indices().begin(),
               edge.adjacent_triangle_indices().end(),
               [this](index_type const f) { return is_boundary_triangle(f); }) == 2u;
}

bool topological_simulated_tetrahedral_mesh_t::is_interior_edge(index_type e) const
{
    edge_t const& edge = edges_[e];
    return std::none_of(
        edge.adjacent_triangle_indices().begin(),
        edge.adjacent_triangle_indices().end(),
        [this](index_type const f) { return is_boundary_triangle(f); });
}

bool topological_simulated_tetrahedral_mesh_t::is_boundary_triangle(index_type f) const
{
    triangle_t const& triangle = triangles_[f];
    return triangle.adjacent_tetrahedron_indices().size() == 1u;
}

bool topological_simulated_tetrahedral_mesh_t::is_interior_triangle(index_type triangle) const
{
    return triangles_[triangle].is_interior_triangle();
}

renderable_topological_simulated_tetrahedral_mesh_t::
    renderable_topological_simulated_tetrahedral_mesh_t(
        common::geometry_t const& geometry,
        build_topology_parameters_t const& params)
    : topological_simulated_tetrahedral_mesh_t(geometry, params)
{
}

void renderable_topological_simulated_tetrahedral_mesh_t::prepare_vertices_for_rendering()
{
    std::size_t constexpr num_vertices_per_triangle = 3u;
    std::size_t constexpr num_attributes_per_vertex = 9u;
    std::size_t const num_boundary_triangles        = boundary_triangle_count();
    std::size_t const vertex_count = num_boundary_triangles * num_vertices_per_triangle;
    std::vector<float> vertex_buffer{};
    vertex_buffer.reserve(vertex_count * num_attributes_per_vertex);

    for (triangle_t const& triangle : triangles())
    {
        if (triangle.is_interior_triangle())
            continue;

        auto const v1 = triangle.v1();
        auto const v2 = triangle.v2();
        auto const v3 = triangle.v3();
        std::array<index_type, 3u> v{v1, v2, v3};

        common::triangle_t triangle_primitive{
            vertices().at(v1).position(),
            vertices().at(v2).position(),
            vertices().at(v3).position()};

        auto const normal = triangle_primitive.normal();

        for (std::size_t i = 0u; i < v.size(); ++i)
        {
            vertex_t const& vertex = vertices().at(v[i]);

            vertex_buffer.push_back(static_cast<float>(vertex.position().x()));
            vertex_buffer.push_back(static_cast<float>(vertex.position().y()));
            vertex_buffer.push_back(static_cast<float>(vertex.position().z()));

            // use triangle normal since face-based
            vertex_buffer.push_back(static_cast<float>(normal.x()));
            vertex_buffer.push_back(static_cast<float>(normal.y()));
            vertex_buffer.push_back(static_cast<float>(normal.z()));

            vertex_buffer.push_back(vertex.color().x());
            vertex_buffer.push_back(vertex.color().y());
            vertex_buffer.push_back(vertex.color().z());
        }
    }

    transfer_vertices_for_rendering(std::move(vertex_buffer));
}

void renderable_topological_simulated_tetrahedral_mesh_t::prepare_indices_for_rendering()
{
    // Can't extract boundary surface if triangle-to-tetrahedra adjacency information
    // does not exist
    if (!topology_parameters().triangle_to_tetrahedra)
        return;

    std::size_t const vertex_count = boundary_triangle_count() * 3u;
    std::vector<std::uint32_t> index_buffer(vertex_count);
    std::iota(index_buffer.begin(), index_buffer.end(), 0u);

    transfer_indices_for_rendering(std::move(index_buffer));
}

tetrahedral_mesh_surface_mesh_adapter_t::tetrahedral_mesh_surface_mesh_adapter_t(
    topological_simulated_tetrahedral_mesh_t const* mesh)
    : mesh_(mesh), index_map_{}, vertices_{}, triangles_{}
{
    extract_boundary_surface();
    extract_surface_normals();
}

std::size_t tetrahedral_mesh_surface_mesh_adapter_t::triangle_count() const
{
    return triangles_.size();
}

std::size_t tetrahedral_mesh_surface_mesh_adapter_t::vertex_count() const
{
    return vertices_.size();
}

common::shared_vertex_surface_mesh_i::vertex_type
tetrahedral_mesh_surface_mesh_adapter_t::vertex(std::size_t vi) const
{
    return vertices_[vi];
}

common::shared_vertex_surface_mesh_i::triangle_type
tetrahedral_mesh_surface_mesh_adapter_t::triangle(std::size_t fi) const
{
    return triangles_[fi];
}

std::vector<index_type> const&
tetrahedral_mesh_surface_mesh_adapter_t::surface_to_tetrahedral_mesh_index_map() const
{
    return index_map_;
}

index_type tetrahedral_mesh_surface_mesh_adapter_t::from_surface_vertex(std::size_t vi) const
{
    return index_map_[vi];
}

void tetrahedral_mesh_surface_mesh_adapter_t::extract_boundary_surface()
{
    index_map_.clear();
    vertices_.clear();
    triangles_.clear();

    std::size_t const tet_mesh_vertex_count = mesh_->vertices().size();

    // maps tetrahedral mesh vertices to surface mesh vertices
    std::vector<std::optional<index_type>> tet_to_surface_index_map(tet_mesh_vertex_count);

    // pre-allocate vertex storage heuristically, and triangle storage exactly
    std::vector<triangle_t const*> boundary_triangles = mesh_->boundary_triangles();
    vertices_.reserve(boundary_triangles.size() / 3u);
    triangles_.reserve(boundary_triangles.size());

    for (auto const& boundary_triangle : boundary_triangles)
    {
        auto const v1 = boundary_triangle->v1();
        auto const v2 = boundary_triangle->v2();
        auto const v3 = boundary_triangle->v3();

        std::array<index_type, 3u> v{v1, v2, v3};

        for (std::size_t j = 0u; j < v.size(); ++j)
        {
            index_type const vi = v[j];

            if (!tet_to_surface_index_map[vi].has_value())
            {
                // update map
                index_type const new_vertex_index = static_cast<index_type>(vertices_.size());
                tet_to_surface_index_map[vi]      = new_vertex_index;

                // create boundary vertex
                vertex_t const& physics_vertex = mesh_->vertices().at(vi);
                vertex_type vertex             = from_physics_vertex(physics_vertex);
                vertices_.push_back(vertex);
            }
        }

        triangle_type triangle{};
        triangle.v1 = tet_to_surface_index_map[v1].value();
        triangle.v2 = tet_to_surface_index_map[v2].value();
        triangle.v3 = tet_to_surface_index_map[v3].value();
        triangles_.push_back(triangle);
    }

    index_map_.resize(vertices_.size());
    for (std::size_t vi = 0u; vi < tet_mesh_vertex_count; ++vi)
    {
        if (!tet_to_surface_index_map[vi].has_value())
            continue;

        index_type const surface_vi = tet_to_surface_index_map[vi].value();
        index_map_[surface_vi]      = static_cast<index_type>(vi);
    }
}

void tetrahedral_mesh_surface_mesh_adapter_t::extract_surface_normals()
{
    std::vector<common::normal_t> normals(vertices_.size(), common::normal_t{0., 0., 0.});

    for (auto const& triangle : triangles_)
    {
        auto const v1 = triangle.v1;
        auto const v2 = triangle.v2;
        auto const v3 = triangle.v3;

        common::triangle_t triangle_primitive{
            common::point_t{vertices_[v1].x, vertices_[v1].y, vertices_[v1].z},
            common::point_t{vertices_[v2].x, vertices_[v2].y, vertices_[v2].z},
            common::point_t{vertices_[v3].x, vertices_[v3].y, vertices_[v3].z}};

        auto const n = triangle_primitive.normal();
        auto const A = triangle_primitive.area();

        std::array<index_type, 3u> const v{v1, v2, v3};
        for (std::size_t j = 0u; j < v.size(); ++j)
        {
            normals[v[j]] += A * n;
        }
    }

    for (std::size_t vi = 0u; vi < vertices_.size(); ++vi)
    {
        vertex_type& vertex      = vertices_[vi];
        common::normal_t const n = normals[vi].normalized();
        vertex.nx                = n.x();
        vertex.ny                = n.y();
        vertex.nz                = n.z();
    }
}

topological_simulated_tetrahedral_mesh_t const*
tetrahedral_mesh_surface_mesh_adapter_t::tetrahedral_mesh() const
{
    return mesh_;
}

tetrahedral_mesh_surface_mesh_adapter_t::vertex_type
tetrahedral_mesh_surface_mesh_adapter_t::from_physics_vertex(vertex_t const& v)
{
    physics::vertex_t const& physics_vertex = v;
    vertex_type vertex{};
    vertex.x  = physics_vertex.position().x();
    vertex.y  = physics_vertex.position().y();
    vertex.z  = physics_vertex.position().z();
    vertex.nx = physics_vertex.normal().x();
    vertex.ny = physics_vertex.normal().y();
    vertex.nz = physics_vertex.normal().z();
    vertex.r  = physics_vertex.color().x();
    vertex.g  = physics_vertex.color().y();
    vertex.b  = physics_vertex.color().z();
    return vertex;
}

} // namespace physics
} // namespace sbs