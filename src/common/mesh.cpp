#include "common/mesh.h"

#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <map>
#include <optional>

namespace sbs {
namespace common {

shared_vertex_mesh_t::shared_vertex_mesh_t(common::geometry_t const& geometry)
{
    auto const num_vertices = geometry.positions.size() / 3u;
    positions_.resize(3u, num_vertices);
    for (std::size_t i = 0u; i < num_vertices; ++i)
    {
        auto const xidx = i * 3u;
        auto const yidx = i * 3u + 1u;
        auto const zidx = i * 3u + 2u;

        auto const x = static_cast<double>(geometry.positions[xidx]);
        auto const y = static_cast<double>(geometry.positions[yidx]);
        auto const z = static_cast<double>(geometry.positions[zidx]);

        positions_.col(i) = Eigen::Vector3d{x, y, z};
    }

    auto const num_colors = geometry.colors.size() / 3u;
    colors_.resize(3u, num_colors);
    for (std::size_t i = 0u; i < num_colors; ++i)
    {
        auto const ridx = i * 3u;
        auto const gidx = i * 3u + 1u;
        auto const bidx = i * 3u + 2u;

        auto const r = static_cast<float>(geometry.colors[ridx]) / 255.f;
        auto const g = static_cast<float>(geometry.colors[gidx]) / 255.f;
        auto const b = static_cast<float>(geometry.colors[bidx]) / 255.f;

        colors_.col(i) = Eigen::Vector3f{r, g, b};
    }

    auto const num_uvs = geometry.uvs.size() / 2u;
    uvs_.resize(2u, num_uvs);
    for (std::size_t i = 0u; i < num_uvs; ++i)
    {
        auto const uidx = i * 2u;
        auto const vidx = i * 2u + 1u;

        auto const u = geometry.uvs[uidx];
        auto const v = geometry.uvs[vidx];

        uvs_.col(i) = Eigen::Vector2f{u, v};
    }

    if (geometry.geometry_type == common::geometry_t::geometry_type_t::triangle)
    {
        auto const num_triangles = geometry.indices.size() / 3u;
        elements_.resize(3u, num_triangles);
        for (std::size_t f = 0u; f < num_triangles; ++f)
        {
            auto const v1idx = f * 3u;
            auto const v2idx = f * 3u + 1u;
            auto const v3idx = f * 3u + 2u;

            auto const v1 = static_cast<std::uint32_t>(geometry.indices[v1idx]);
            auto const v2 = static_cast<std::uint32_t>(geometry.indices[v2idx]);
            auto const v3 = static_cast<std::uint32_t>(geometry.indices[v3idx]);

            elements_.col(f) = triangle_type{v1, v2, v3};
        }
        boundary_vertices_ = positions_;
        boundary_faces_    = elements_;
        boundary_colors_   = colors_;
        boundary_uvs_      = uvs_;
    }

    if (geometry.geometry_type == common::geometry_t::geometry_type_t::tetrahedron)
    {
        auto const num_tets = geometry.indices.size() / 4u;
        elements_.resize(4u, num_tets);
        for (std::size_t e = 0u; e < num_tets; ++e)
        {
            auto const v1idx = e * 4u;
            auto const v2idx = e * 4u + 1u;
            auto const v3idx = e * 4u + 2u;
            auto const v4idx = e * 4u + 3u;

            auto const v1 = static_cast<std::uint32_t>(geometry.indices[v1idx]);
            auto const v2 = static_cast<std::uint32_t>(geometry.indices[v2idx]);
            auto const v3 = static_cast<std::uint32_t>(geometry.indices[v3idx]);
            auto const v4 = static_cast<std::uint32_t>(geometry.indices[v4idx]);

            elements_.col(e) = tetrahedron_type{v1, v2, v3, v4};
        }
    }

    masses_.resize(positions_.cols());
    masses_.setOnes();
    velocities_.resize(3u, positions_.cols());
    velocities_.setZero();
    forces_.resize(3u, positions_.cols());
    forces_.setZero();
}

shared_vertex_mesh_t::shared_vertex_mesh_t(
    positions_type const& P,
    tetrahedra_type const& T)
    : positions_(P), elements_(T)
{
    masses_.resize(P.cols());
    velocities_.resizeLike(P);
    forces_.resizeLike(P);
}

shared_vertex_mesh_t::positions_type const& shared_vertex_mesh_t::positions() const
{
    return positions_;
}

shared_vertex_mesh_t::positions_type& shared_vertex_mesh_t::positions()
{
    return positions_;
}

shared_vertex_mesh_t::elements_type const& shared_vertex_mesh_t::elements() const
{
    return elements_;
}

shared_vertex_mesh_t::elements_type& shared_vertex_mesh_t::elements()
{
    return elements_;
}

shared_vertex_mesh_t::masses_type const& shared_vertex_mesh_t::masses() const
{
    return masses_;
}

shared_vertex_mesh_t::masses_type& shared_vertex_mesh_t::masses()
{
    return masses_;
}

shared_vertex_mesh_t::velocities_type const& shared_vertex_mesh_t::velocities() const
{
    return velocities_;
}

shared_vertex_mesh_t::velocities_type& shared_vertex_mesh_t::velocities()
{
    return velocities_;
}

shared_vertex_mesh_t::forces_type const& shared_vertex_mesh_t::forces() const
{
    return forces_;
}

shared_vertex_mesh_t::forces_type& shared_vertex_mesh_t::forces()
{
    return forces_;
}

shared_vertex_mesh_t::positions_type const& shared_vertex_mesh_t::vertices() const
{
    return boundary_vertices_;
}

shared_vertex_mesh_t::positions_type& shared_vertex_mesh_t::vertices()
{
    return boundary_vertices_;
}

shared_vertex_mesh_t::uv_coordinates_type const& shared_vertex_mesh_t::uvs() const
{
    return boundary_uvs_;
}
shared_vertex_mesh_t::uv_coordinates_type& shared_vertex_mesh_t::uvs()
{
    return boundary_uvs_;
}

shared_vertex_mesh_t::normals_type const& shared_vertex_mesh_t::normals() const
{
    return normals_;
}
shared_vertex_mesh_t::normals_type& shared_vertex_mesh_t::normals()
{
    return normals_;
}

shared_vertex_mesh_t::colors_type const& shared_vertex_mesh_t::colors() const
{
    return boundary_colors_;
}
shared_vertex_mesh_t::colors_type& shared_vertex_mesh_t::colors()
{
    return boundary_colors_;
}

void shared_vertex_mesh_t::set_color(Eigen::Vector3f const rgb)
{
    colors_.resizeLike(positions_);
    colors_.row(0u).setConstant(rgb(0u));
    colors_.row(1u).setConstant(rgb(1u));
    colors_.row(2u).setConstant(rgb(2u));

    boundary_colors_.row(0u).setConstant(rgb(0u));
    boundary_colors_.row(1u).setConstant(rgb(1u));
    boundary_colors_.row(2u).setConstant(rgb(2u));
}

shared_vertex_mesh_t::triangles_type const& shared_vertex_mesh_t::faces() const
{
    return boundary_faces_;
}

shared_vertex_mesh_t::triangles_type& shared_vertex_mesh_t::faces()
{
    return boundary_faces_;
}

void shared_vertex_mesh_t::extract_boundary_surface_mesh()
{
    bool const is_triangle_mesh = elements_.rows() == 3u;
    bool const is_tet_mesh      = elements_.rows() == 4u;

    /**
     * This is a triangle mesh
     */
    if (is_triangle_mesh)
    {
        boundary_vertices_ = positions_;
        boundary_faces_    = elements_;
        boundary_colors_   = colors_;
        boundary_uvs_      = uvs_;
        return;
    }

    if (is_tet_mesh)
    {
        triangles_type triangles{};
        triangles.resize(3u, elements_.cols() * 4u);

        /**
         * Extract all triangles from tet mesh (all 4 faces of all tets)
         */
        index_type const num_tetrahedra = static_cast<index_type>(elements_.cols());
        for (index_type e = 0u; e < num_tetrahedra; ++e)
        {
            tetrahedron_type const tetrahedron = elements_.col(e);
            auto const v1                      = tetrahedron(0u);
            auto const v2                      = tetrahedron(1u);
            auto const v3                      = tetrahedron(2u);
            auto const v4                      = tetrahedron(3u);

            index_type const f1 = e * 4u;
            index_type const f2 = e * 4u + 1u;
            index_type const f3 = e * 4u + 2u;
            index_type const f4 = e * 4u + 3u;

            triangles(0u, f1) = v1;
            triangles(1u, f1) = v3;
            triangles(2u, f1) = v2;

            triangles(0u, f2) = v1;
            triangles(1u, f2) = v2;
            triangles(2u, f2) = v4;

            triangles(0u, f3) = v2;
            triangles(1u, f3) = v3;
            triangles(2u, f3) = v4;

            triangles(0u, f4) = v1;
            triangles(1u, f4) = v4;
            triangles(2u, f4) = v3;
        }

        /**
         * Reorder triangle indices in the same order for all triangles
         */
        std::vector<std::array<index_type, 3u>> triangle_copies{};
        triangle_copies.reserve(triangles.cols());

        index_type const num_triangles = static_cast<index_type>(triangles.cols());
        for (std::size_t f = 0u; f < num_triangles; ++f)
        {
            std::array<index_type, 3u> triangle{
                triangles(0u, f),
                triangles(1u, f),
                triangles(2u, f)};

            std::sort(triangle.begin(), triangle.end());
            triangle_copies.push_back(triangle);
        }

        /**
         * Count number of times each triangle is present in the tet mesh.
         */
        std::map<std::array<index_type, 3u>, std::uint32_t> triangle_occurrences{};
        for (auto const& triangle : triangle_copies)
        {
            if (triangle_occurrences.find(triangle) == triangle_occurrences.end())
            {
                triangle_occurrences[triangle] = 1u;
            }
            else
            {
                ++triangle_occurrences[triangle];
            }
        }

        std::vector<std::array<index_type, 3u>> boundary_triangles{};
        boundary_triangles.reserve(triangle_copies.size());

        std::vector<std::optional<index_type>> index_map(
            positions_.cols(),
            std::optional<index_type>{});

        index_type boundary_vertices_size = 0u;
        for (std::size_t f = 0u; f < triangle_copies.size(); ++f)
        {
            /**
             * Boundary triangles should only be present once in the tet mesh, since
             * there is no adjacent face to it.
             */
            if (triangle_occurrences[triangle_copies[f]] != 1u)
                continue;

            auto const v1 = triangles(0u, f);
            auto const v2 = triangles(1u, f);
            auto const v3 = triangles(2u, f);

            if (!index_map[v1].has_value())
            {
                index_map[v1] = boundary_vertices_size++;
            }
            if (!index_map[v2].has_value())
            {
                index_map[v2] = boundary_vertices_size++;
            }
            if (!index_map[v3].has_value())
            {
                index_map[v3] = boundary_vertices_size++;
            }

            boundary_triangles.push_back(
                {index_map[v1].value(), index_map[v2].value(), index_map[v3].value()});
        }

        bool const has_uvs = uvs_.cols() != 0u;

        boundary_vertices_.resize(3u, boundary_vertices_size);
        boundary_colors_.resize(3u, boundary_vertices_size);
        boundary_uvs_.resize(2u, boundary_vertices_size);

        for (std::size_t i = 0u; i < index_map.size(); ++i)
        {
            if (index_map[i].has_value())
            {
                boundary_vertices_.col(index_map[i].value()) = positions_.col(i);
                boundary_colors_.col(index_map[i].value())   = colors_.col(i);

                if (has_uvs)
                    boundary_uvs_.col(index_map[i].value()) = uvs_.col(i);
            }
        }

        boundary_faces_.resize(3u, boundary_triangles.size());
        for (std::size_t f = 0u; f < boundary_triangles.size(); ++f)
        {
            boundary_faces_(0u, f) = boundary_triangles[f][0u];
            boundary_faces_(1u, f) = boundary_triangles[f][1u];
            boundary_faces_(2u, f) = boundary_triangles[f][2u];
        }
    }
}

void shared_vertex_mesh_t::extract_boundary_normals()
{
    normals_.resizeLike(boundary_vertices_);
    normals_.setZero();

    index_type const num_boundary_faces = static_cast<index_type>(boundary_faces_.cols());
    for (std::size_t f = 0u; f < num_boundary_faces; ++f)
    {
        triangle_type const triangle = boundary_faces_.col(f);

        auto const v1 = triangle(0u);
        auto const v2 = triangle(1u);
        auto const v3 = triangle(2u);

        Eigen::Vector3d const p1 = boundary_vertices_.col(v1);
        Eigen::Vector3d const p2 = boundary_vertices_.col(v2);
        Eigen::Vector3d const p3 = boundary_vertices_.col(v3);

        Eigen::Vector3d const p21 = p2 - p1;
        Eigen::Vector3d const p31 = p3 - p1;

        Eigen::Vector3d const area_weighted_normal = 0.5 * p21.cross(p31);
        normals_.col(v1) += area_weighted_normal;
        normals_.col(v2) += area_weighted_normal;
        normals_.col(v3) += area_weighted_normal;
    }

    normals_.colwise().normalize();
}

void shared_vertex_mesh_t::rescale(Eigen::Vector3d const& boxmin, Eigen::Vector3d const& boxmax)
{
    Eigen::Vector3d const min = positions_.rowwise().minCoeff();
    Eigen::Vector3d const max = positions_.rowwise().maxCoeff();

    double const dx = max.x() - min.x();
    double const dy = max.y() - min.y();
    double const dz = max.z() - min.z();

    double constexpr eps           = 1e-8;
    bool const dx_division_by_zero = std::abs(dx) < eps;
    bool const dy_division_by_zero = std::abs(dy) < eps;
    bool const dz_division_by_zero = std::abs(dz) < eps;

    index_type const num_positions = static_cast<index_type>(positions_.cols());
    for (std::size_t i = 0u; i < num_positions; ++i)
    {
        Eigen::Vector3d const p = positions_.col(i);

        auto const x = dx_division_by_zero ?
                           p.x() :
                           boxmin.x() + (boxmax.x() - boxmin.x()) * (p.x() - min.x()) / dx;
        auto const y = dy_division_by_zero ?
                           p.y() :
                           boxmin.y() + (boxmax.y() - boxmin.y()) * (p.y() - min.y()) / dy;
        auto const z = dz_division_by_zero ?
                           p.z() :
                           boxmin.z() + (boxmax.z() - boxmin.z()) * (p.z() - min.z()) / dz;

        positions_(0u, i) = x;
        positions_(1u, i) = y;
        positions_(2u, i) = z;
    }

    index_type const num_boundary_vertices = static_cast<index_type>(boundary_vertices_.cols());
    for (std::size_t i = 0u; i < num_boundary_vertices; ++i)
    {
        Eigen::Vector3d const p = boundary_vertices_.col(i);

        auto const x = dx_division_by_zero ?
                           p.x() :
                           boxmin.x() + (boxmax.x() - boxmin.x()) * (p.x() - min.x()) / dx;
        auto const y = dy_division_by_zero ?
                           p.y() :
                           boxmin.y() + (boxmax.y() - boxmin.y()) * (p.y() - min.y()) / dy;
        auto const z = dz_division_by_zero ?
                           p.z() :
                           boxmin.z() + (boxmax.z() - boxmin.z()) * (p.z() - min.z()) / dz;

        boundary_vertices_(0u, i) = x;
        boundary_vertices_(1u, i) = y;
        boundary_vertices_(2u, i) = z;
    }
}

std::vector<std::pair<std::uint32_t, std::uint32_t>> edges(shared_vertex_mesh_t const& mesh)
{
    using edge_type = std::pair<std::uint32_t, std::uint32_t>;

    std::vector<edge_type> edges_with_duplicates{};

    std::size_t const num_elements              = static_cast<std::size_t>(mesh.elements().cols());
    std::size_t constexpr num_edges_per_element = 6u;
    edges_with_duplicates.reserve(num_elements * num_edges_per_element);

    using element_type = Eigen::Matrix<std::uint32_t, Eigen::Dynamic, 1>;
    for (std::size_t e = 0u; e < num_elements; ++e)
    {
        element_type const element = mesh.elements().col(e);
        /**
         * Triangle
         */
        if (element.rows() == 3u)
        {
            auto const v1 = element(0u);
            auto const v2 = element(1u);
            auto const v3 = element(2u);

            /**
             * Add index pairs in sorted order
             */
            edges_with_duplicates.push_back(
                v1 < v2 ? std::make_pair(v1, v2) : std::make_pair(v2, v1));
            edges_with_duplicates.push_back(
                v2 < v3 ? std::make_pair(v2, v3) : std::make_pair(v3, v2));
            edges_with_duplicates.push_back(
                v3 < v1 ? std::make_pair(v3, v1) : std::make_pair(v1, v3));
        }
        /**
         * Tetrahedron
         */
        if (element.rows() == 4u)
        {
            auto const v1 = element(0u);
            auto const v2 = element(1u);
            auto const v3 = element(2u);
            auto const v4 = element(3u);

            /**
             * Add index pairs in sorted order
             */
            edges_with_duplicates.push_back(
                v1 < v2 ? std::make_pair(v1, v2) : std::make_pair(v2, v1));
            edges_with_duplicates.push_back(
                v2 < v3 ? std::make_pair(v2, v3) : std::make_pair(v3, v2));
            edges_with_duplicates.push_back(
                v3 < v1 ? std::make_pair(v3, v1) : std::make_pair(v1, v3));

            edges_with_duplicates.push_back(
                v1 < v4 ? std::make_pair(v1, v4) : std::make_pair(v4, v1));
            edges_with_duplicates.push_back(
                v2 < v4 ? std::make_pair(v2, v4) : std::make_pair(v4, v2));
            edges_with_duplicates.push_back(
                v3 < v4 ? std::make_pair(v3, v4) : std::make_pair(v4, v3));
        }
    }

    std::sort(
        edges_with_duplicates.begin(),
        edges_with_duplicates.end(),
        [](edge_type const& e1, edge_type const& e2) {
            if (e1.first < e2.first)
                return true;

            if (e1.first > e2.first)
                return false;

            return e1.second < e2.second;
        });

    auto it = std::unique(edges_with_duplicates.begin(), edges_with_duplicates.end());
    return std::vector<edge_type>(edges_with_duplicates.begin(), it);
}

} // namespace common
} // namespace sbs