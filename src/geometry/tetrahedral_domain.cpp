#include "sbs/geometry/tetrahedral_domain.h"

#include <unordered_map>

namespace sbs {
namespace geometry {

tetrahedral_domain_t::tetrahedral_domain_t(
    std::vector<Eigen::Vector3d> const& points,
    std::vector<index_type> const& indices,
    scalar_type const query_error)
    : positions_(points), mesh_(), tet_maps_()
{
    mesh_.reserve_vertices(points.size());
    for (auto const& p : points)
        mesh_.add_vertex();

    for (auto i = 0u; i < indices.size(); i += 4u)
    {
        index_type const v1 = indices[i];
        index_type const v2 = indices[i + 1u];
        index_type const v3 = indices[i + 2u];
        index_type const v4 = indices[i + 3u];

        mesh_.add_tetrahedron(topology::tetrahedron_t{v1, v2, v3, v4});

        Eigen::Vector3d const& p1 = points[v1];
        Eigen::Vector3d const& p2 = points[v2];
        Eigen::Vector3d const& p3 = points[v3];
        Eigen::Vector3d const& p4 = points[v4];

        tet_maps_.push_back(math::tetrahedron_barycentric_mapping_t(p1, p2, p3, p4));
    }

    in_tetrahedron_query_ = in_tetrahedron_query_t(&mesh_, &positions_, &tet_maps_, query_error);
}

topology::tetrahedron_t const& tetrahedral_domain_t::tetrahedron(index_type ti) const
{
    return mesh_.tetrahedron(ti);
}

Eigen::Vector3d const& tetrahedral_domain_t::position(index_type i) const
{
    return positions_[i];
}

topology::tetrahedron_set_t const& tetrahedral_domain_t::topology() const
{
    return mesh_;
}

index_type tetrahedral_domain_t::in_tetrahedron(Eigen::Vector3d const& X) const
{
    return in_tetrahedron_query_.in_tetrahedron(X);
}

math::tetrahedron_barycentric_mapping_t const&
tetrahedral_domain_t::barycentric_map(index_type ti) const
{
    return tet_maps_[ti];
}

tetrahedral_domain_t::in_tetrahedron_query_t::in_tetrahedron_query_t()
    : base_type(0u),
      topology_(nullptr),
      mesh_nodes_(nullptr),
      tet_maps_(nullptr),
      tetrahedron_centers_()
{
}

tetrahedral_domain_t::in_tetrahedron_query_t::in_tetrahedron_query_t(
    topology::tetrahedron_set_t const* topology,
    std::vector<Eigen::Vector3d> const* mesh_nodes,
    std::vector<math::tetrahedron_barycentric_mapping_t> const* tet_maps,
    scalar_type tolerance)
    : base_type(topology->tetrahedron_count()),
      topology_(topology),
      mesh_nodes_(mesh_nodes),
      tet_maps_(tet_maps),
      tetrahedron_centers_(),
      tolerance_(tolerance)
{
    tetrahedron_centers_.reserve(topology_->tetrahedron_count());
    for (auto const& t : topology_->tetrahedra())
    {
        Eigen::Vector3d const& p1 = (*mesh_nodes)[t.v1()];
        Eigen::Vector3d const& p2 = (*mesh_nodes)[t.v2()];
        Eigen::Vector3d const& p3 = (*mesh_nodes)[t.v3()];
        Eigen::Vector3d const& p4 = (*mesh_nodes)[t.v4()];

        Eigen::Vector3d const center = 0.25 * (p1 + p2 + p3 + p4);
        tetrahedron_centers_.push_back(center);
    }
    this->construct();
}

index_type
tetrahedral_domain_t::in_tetrahedron_query_t::in_tetrahedron(Eigen::Vector3d const& p) const
{
    auto const intersects = [this, p](unsigned int node_idx, unsigned int depth) -> bool {
        Discregrid::BoundingSphere const& s = this->hull(node_idx);
        return s.contains(p);
    };

    auto const is_point_in_tetrahedron =
        [this](Eigen::Vector3d const& point, topology::tetrahedron_t const& t) {
            auto const& face_copies = t.faces_copy();
            std::array<bool, 4u> is_inside{false, false, false, false};
            for (std::uint8_t i = 0u; i < 4u; ++i)
            {
                auto const& f             = face_copies[i];
                Eigen::Vector3d const& p1 = (*mesh_nodes_)[f.v1()];
                Eigen::Vector3d const& p2 = (*mesh_nodes_)[f.v2()];
                Eigen::Vector3d const& p3 = (*mesh_nodes_)[f.v3()];

                Eigen::Vector3d const n           = (p2 - p1).cross(p3 - p1).normalized();
                scalar_type const signed_distance = (point - p1).dot(n);
                bool const is_outside             = signed_distance > sbs::eps();
                is_inside[i]                      = !is_outside;
            }
            bool const is_p_in_t = is_inside[0] && is_inside[1] && is_inside[2] && is_inside[3];
            return is_p_in_t;
        };

    index_type parent_ti = std::numeric_limits<index_type>::max();
    bool found{false};
    auto const get_ti = [this, p, &parent_ti, is_point_in_tetrahedron, &found](
                            unsigned int node_idx,
                            unsigned int depth) {
        if (found)
            return;

        base_type::Node const& node = this->node(node_idx);
        if (!node.isLeaf())
            return;

        for (auto j = node.begin; j < node.begin + node.n; ++j)
        {
            index_type const ti              = static_cast<index_type>(m_lst[j]);
            topology::tetrahedron_t const& t = topology_->tetrahedron(ti);

            Eigen::Vector3d const& p1 = (*mesh_nodes_)[t.v1()];
            Eigen::Vector3d const& p2 = (*mesh_nodes_)[t.v2()];
            Eigen::Vector3d const& p3 = (*mesh_nodes_)[t.v3()];
            Eigen::Vector3d const& p4 = (*mesh_nodes_)[t.v4()];

            bool const is_point_contained = is_point_in_tetrahedron(p, t);

            if (is_point_contained)
            {
                // NOTE:
                // Use this variable for debugging purposes.
                // It is possible that a point lies exactly in 2 tetrahedra or more.
                // For example, a point lying on a shared face of 2 tetrahedra will
                // belong in both tets. In this implementation, for the moment,
                // we will choose to associate a point with the first tetrahedra
                // that we find.
                bool const is_point_shared_between_multiple_tetrahedra =
                    (parent_ti != std::numeric_limits<index_type>::max());
                parent_ti = ti;
                found     = true;
                break;
            }
        }
    };

    traverseBreadthFirst(intersects, get_ti);

    return parent_ti;
}

Eigen::Vector3d tetrahedral_domain_t::in_tetrahedron_query_t::entityPosition(unsigned int i) const
{
    return tetrahedron_centers_[i];
}

void tetrahedral_domain_t::in_tetrahedron_query_t::computeHull(
    unsigned int b,
    unsigned int n,
    Discregrid::BoundingSphere& hull) const
{
    std::vector<Eigen::Vector3d> vertices_of_sphere{};
    vertices_of_sphere.reserve(n * 4u);
    for (unsigned int i = b; i < n + b; ++i)
    {
        index_type const ti              = static_cast<index_type>(m_lst[i]);
        topology::tetrahedron_t const& t = topology_->tetrahedron(ti);
        for (index_type const vi : t.vertex_indices())
        {
            Eigen::Vector3d const& pi = (*mesh_nodes_)[vi];
            vertices_of_sphere.push_back(pi);
        }
    }

    Discregrid::BoundingSphere const s(vertices_of_sphere);

    hull.x() = s.x();
    hull.r() = s.r() + tolerance_;
}

std::pair<std::vector<Eigen::Vector3d>, std::vector<index_type>>
boundary_surface(tetrahedral_domain_t const& domain)
{
    auto const& topology = domain.topology();
    std::vector<Eigen::Vector3d> tet_points{};
    tet_points.reserve(topology.vertex_count());
    for (auto i = 0u; i < topology.vertex_count(); ++i)
    {
        tet_points.push_back(domain.position(i));
    }

    std::vector<topology::triangle_t> const boundary_triangles =
        topology.oriented_boundary_triangles();

    std::vector<Eigen::Vector3d> triangle_points{};
    triangle_points.reserve(tet_points.size()); // overallocate

    std::vector<index_type> triangle_indices{};
    triangle_indices.reserve(boundary_triangles.size()); // heuristically pre-allocate

    std::unordered_map<index_type, index_type> tet_to_surface_vertex_map{};
    for (auto const& f : boundary_triangles)
    {
        if (tet_to_surface_vertex_map.find(f.v1()) == tet_to_surface_vertex_map.end())
        {
            auto const index                  = tet_to_surface_vertex_map.size();
            tet_to_surface_vertex_map[f.v1()] = static_cast<index_type>(index);
            triangle_points.push_back(tet_points[f.v1()]);
        }
        if (tet_to_surface_vertex_map.find(f.v2()) == tet_to_surface_vertex_map.end())
        {
            auto const index                  = tet_to_surface_vertex_map.size();
            tet_to_surface_vertex_map[f.v2()] = static_cast<index_type>(index);
            triangle_points.push_back(tet_points[f.v2()]);
        }
        if (tet_to_surface_vertex_map.find(f.v3()) == tet_to_surface_vertex_map.end())
        {
            auto const index                  = tet_to_surface_vertex_map.size();
            tet_to_surface_vertex_map[f.v3()] = static_cast<index_type>(index);
            triangle_points.push_back(tet_points[f.v3()]);
        }

        auto const v1 = tet_to_surface_vertex_map[f.v1()];
        auto const v2 = tet_to_surface_vertex_map[f.v2()];
        auto const v3 = tet_to_surface_vertex_map[f.v3()];
        triangle_indices.push_back(v1);
        triangle_indices.push_back(v2);
        triangle_indices.push_back(v3);
    }

    return std::make_pair(triangle_points, triangle_indices);
}

} // namespace geometry
} // namespace sbs