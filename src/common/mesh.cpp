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
        boundary_faces_.resize(3u, num_triangles);
        for (std::size_t f = 0u; f < num_triangles; ++f)
        {
            auto const v1idx = f * 3u;
            auto const v2idx = f * 3u + 1u;
            auto const v3idx = f * 3u + 2u;

            auto const v1 = static_cast<std::uint32_t>(geometry.indices[v1idx]);
            auto const v2 = static_cast<std::uint32_t>(geometry.indices[v2idx]);
            auto const v3 = static_cast<std::uint32_t>(geometry.indices[v3idx]);

            boundary_faces_.col(f) = triangle_type{v1, v2, v3};
        }
        boundary_vertices_ = positions_;
        boundary_colors_   = colors_;
    }

    if (geometry.geometry_type == common::geometry_t::geometry_type_t::tetrahedron)
    {
        auto const num_tets = geometry.indices.size() / 4u;
        tetrahedra_.resize(4u, num_tets);
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

            tetrahedra_.col(e) = tetrahedron_type{v1, v2, v3, v4};
        }
        this->extract_boundary_surface_mesh();
    }

    this->extract_boundary_normals();
}

shared_vertex_mesh_t::positions_type const& shared_vertex_mesh_t::positions() const
{
    return positions_;
}

shared_vertex_mesh_t::positions_type& shared_vertex_mesh_t::positions()
{
    return positions_;
}

shared_vertex_mesh_t::tetrahedra_type const& shared_vertex_mesh_t::elements() const
{
    return tetrahedra_;
}

shared_vertex_mesh_t::tetrahedra_type& shared_vertex_mesh_t::elements()
{
    return tetrahedra_;
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

shared_vertex_mesh_t::positions_type const& shared_vertex_mesh_t::boundary_vertices() const
{
    return boundary_vertices_;
}

shared_vertex_mesh_t::positions_type& shared_vertex_mesh_t::boundary_vertices()
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
    triangles_type triangles{};
    triangles.resize(3u, tetrahedra_.cols() * 4u);

    /**
     * Extract all triangles from tet mesh (all 4 faces of all tets)
     */
    for (std::size_t e = 0u; e < tetrahedra_.cols(); ++e)
    {
        tetrahedron_type const tetrahedron = tetrahedra_.col(e);
        auto const v1                      = tetrahedron(0u);
        auto const v2                      = tetrahedron(1u);
        auto const v3                      = tetrahedron(2u);
        auto const v4                      = tetrahedron(3u);

        std::size_t const f1 = e * 4u;
        std::size_t const f2 = e * 4u + 1u;
        std::size_t const f3 = e * 4u + 2u;
        std::size_t const f4 = e * 4u + 3u;

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
    std::vector<std::array<std::uint32_t, 3u>> triangle_copies{};
    triangle_copies.reserve(triangles.cols());

    for (std::size_t f = 0u; f < triangles.cols(); ++f)
    {
        std::array<std::uint32_t, 3u> triangle{
            triangles(0u, f),
            triangles(1u, f),
            triangles(2u, f)};

        std::sort(triangle.begin(), triangle.end());
        triangle_copies.push_back(triangle);
    }

    /**
     * Count number of times each triangle is present in the tet mesh.
     */
    std::map<std::array<std::uint32_t, 3u>, std::uint32_t> triangle_occurrences{};
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

    std::vector<std::array<std::uint32_t, 3u>> boundary_triangles{};
    boundary_triangles.reserve(triangle_copies.size());

    std::vector<std::optional<std::uint32_t>> index_map(
        positions_.cols(),
        std::optional<std::uint32_t>{});

    std::size_t boundary_vertices_size = 0u;
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

void shared_vertex_mesh_t::extract_boundary_normals()
{
    normals_.resizeLike(boundary_vertices_);
    normals_.setZero();

    for (std::size_t f = 0u; f < boundary_faces_.cols(); ++f)
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

    for (std::size_t i = 0u; i < positions_.cols(); ++i)
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

    for (std::size_t i = 0u; i < boundary_vertices_.cols(); ++i)
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

} // namespace common
} // namespace sbs