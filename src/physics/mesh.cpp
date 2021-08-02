#include "sbs/physics/mesh.h"

namespace sbs {
namespace physics {

vertex_t::vertex_t()
    : p_{0., 0., 0.},
      v_{0., 0., 0.},
      f_{0., 0., 0.},
      m_{0.},
      n_{0., 0., 0.},
      c_{0.f, 0.f, 0.f},
      fixed_{false}
{
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

bool vertex_t::fixed() const
{
    return fixed_;
}

bool& vertex_t::fixed()
{
    return fixed_;
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

tetrahedron_t& tetrahedron_t::operator=(tetrahedron_t const& other)
{
    v_   = other.v_;
    rho_ = other.rho_;
    if (static_cast<bool>(other.edges_))
    {
        edges_ = std::make_unique<std::array<index_type, 6u>>(*other.edges_);
    }
    if (static_cast<bool>(other.faces_))
    {
        faces_ = std::make_unique<std::array<index_type, 4u>>(*other.faces_);
    }
    return *this;
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

std::array<index_type, 4u> const& tetrahedron_t::vertices() const
{
    return v_;
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

std::array<index_type, 2u> const& edge_t::vertices() const
{
    return v_;
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

triangle_t& triangle_t::operator=(triangle_t const& other)
{
    v_             = other.v_;
    adjacent_tets_ = other.adjacent_tets_;
    if (static_cast<bool>(other.edges_))
    {
        edges_ = std::make_unique<std::array<index_type, 3u>>(*other.edges_);
    }
    return *this;
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

std::array<index_type, 3u> const& triangle_t::vertices() const
{
    return v_;
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

} // namespace physics
} // namespace sbs