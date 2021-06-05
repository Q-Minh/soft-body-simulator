#include "common/mesh.h"

#include <Eigen/Geometry>
#include <algorithm>
#include <array>
#include <map>
#include <numeric>
#include <optional>

namespace sbs {
namespace common {

shared_vertex_mesh_t::shared_vertex_mesh_t(common::geometry_t const& geometry)
{
    if (geometry.geometry_type != common::geometry_t::geometry_type_t::tetrahedron)
        return;

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

    masses_.resize(positions_.cols());
    masses_.setOnes();
    velocities_.resize(3u, positions_.cols());
    velocities_.setZero();
    forces_.resize(3u, positions_.cols());
    forces_.setZero();
}

shared_vertex_mesh_t::shared_vertex_mesh_t(positions_type const& P, tetrahedra_type const& T)
    : positions_(P), elements_(T)
{
    masses_.resize(P.cols());
    velocities_.resizeLike(P);
    forces_.resizeLike(P);
}

shared_vertex_mesh_t::shared_vertex_mesh_t(
    positions_type const& P,
    tetrahedra_type const& T,
    masses_type const& M,
    velocities_type const& V,
    forces_type const& F)
    : positions_(P), elements_(T), masses_(M), velocities_(V), forces_(F)
{
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

shared_vertex_surface_mesh_t
shared_vertex_mesh_t::boundary_surface_mesh(Eigen::Vector3f const& color) const
{
    shared_vertex_surface_mesh_t surface_mesh{};
    surface_mesh.triangles().resize(3u, elements_.cols() * 4u);

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

        surface_mesh.triangles()(0u, f1) = v1;
        surface_mesh.triangles()(1u, f1) = v3;
        surface_mesh.triangles()(2u, f1) = v2;

        surface_mesh.triangles()(0u, f2) = v1;
        surface_mesh.triangles()(1u, f2) = v2;
        surface_mesh.triangles()(2u, f2) = v4;

        surface_mesh.triangles()(0u, f3) = v2;
        surface_mesh.triangles()(1u, f3) = v3;
        surface_mesh.triangles()(2u, f3) = v4;

        surface_mesh.triangles()(0u, f4) = v1;
        surface_mesh.triangles()(1u, f4) = v4;
        surface_mesh.triangles()(2u, f4) = v3;
    }

    /**
     * Reorder triangle indices in the same order for all triangles
     */
    std::vector<std::array<index_type, 3u>> triangle_copies{};
    triangle_copies.reserve(surface_mesh.triangles().cols());

    index_type const num_triangles = static_cast<index_type>(surface_mesh.triangles().cols());
    for (std::size_t f = 0u; f < num_triangles; ++f)
    {
        std::array<index_type, 3u> triangle{
            surface_mesh.triangles()(0u, f),
            surface_mesh.triangles()(1u, f),
            surface_mesh.triangles()(2u, f)};

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

    std::size_t num_boundary_triangles  = 0u;
    std::uint32_t num_boundary_vertices = 0u;
    /**
     * Transfer boundary triangles to surface mesh while building an index map
     * which maps indices to the vertices of the tet mesh to indices to the
     * vertices of the triangle mesh, since we rebuild our vertices buffer
     * to discard all interior tet vertices. We don't want to need to transfer
     * useless vertices to the GPU when rendering.
     */
    for (std::size_t f = 0u; f < triangle_copies.size(); ++f)
    {
        /**
         * Boundary triangles should only be present once in the tet mesh, since
         * there is no adjacent face to it.
         */
        if (triangle_occurrences[triangle_copies[f]] != 1u)
            continue;

        auto const v1 = surface_mesh.triangles()(0u, f);
        auto const v2 = surface_mesh.triangles()(1u, f);
        auto const v3 = surface_mesh.triangles()(2u, f);

        if (!index_map[v1].has_value())
        {
            index_map[v1] = num_boundary_vertices++;
        }
        if (!index_map[v2].has_value())
        {
            index_map[v2] = num_boundary_vertices++;
        }
        if (!index_map[v3].has_value())
        {
            index_map[v3] = num_boundary_vertices++;
        }

        surface_mesh.triangles().col(num_boundary_triangles++) =
            shared_vertex_surface_mesh_t::triangle_type{
                index_map[v1].value(),
                index_map[v2].value(),
                index_map[v3].value()};
    }
    surface_mesh.triangles().conservativeResize(3u, num_boundary_triangles);

    /**
     * Transfer only boundary vertices to the surface mesh
     */
    surface_mesh.vertices().resize(3u, num_boundary_vertices);
    surface_mesh.index_map().resize(num_boundary_vertices);
    for (std::size_t i = 0u; i < index_map.size(); ++i)
    {
        if (!index_map[i].has_value())
            continue;

        auto const surface_mesh_vertex_index                   = index_map[i].value();
        surface_mesh.index_map()[surface_mesh_vertex_index]    = static_cast<index_type>(i);
        surface_mesh.vertices().col(surface_mesh_vertex_index) = positions_.col(i);
    }

    surface_mesh.extract_normals();
    surface_mesh.set_color(color);
    return surface_mesh;
}

shared_vertex_surface_mesh_t shared_vertex_mesh_t::facets(Eigen::Vector3f const& color) const
{
    shared_vertex_surface_mesh_t triangle_mesh{};
    triangle_mesh.vertices() = positions_;

    /**
     * All vertex indices from the triangle mesh and the original tetrahedral mesh
     * correspond one to one.
     */
    triangle_mesh.index_map().resize(static_cast<std::size_t>(positions_.cols()));
    std::iota(triangle_mesh.index_map().begin(), triangle_mesh.index_map().end(), 0u);

    auto& triangles = triangle_mesh.triangles();
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

    triangle_mesh.extract_normals();
    triangle_mesh.set_color(color);
    return triangle_mesh;
}

shared_vertex_surface_mesh_t::shared_vertex_surface_mesh_t(common::geometry_t const& geometry)
{
    if (geometry.geometry_type != geometry_t::geometry_type_t::triangle)
        return;

    auto const num_vertices = geometry.positions.size() / 3u;
    vertices_.resize(3u, num_vertices);
    for (std::size_t i = 0u; i < num_vertices; ++i)
    {
        auto const xidx = i * 3u;
        auto const yidx = i * 3u + 1u;
        auto const zidx = i * 3u + 2u;

        auto const x = static_cast<double>(geometry.positions[xidx]);
        auto const y = static_cast<double>(geometry.positions[yidx]);
        auto const z = static_cast<double>(geometry.positions[zidx]);

        vertices_.col(i) = Eigen::Vector3d{x, y, z};
    }

    auto const num_normals = geometry.normals.size() / 3u;
    normals_.resize(3u, num_normals);
    for (std::size_t i = 0u; i < num_normals; ++i)
    {
        auto const nxidx = i * 3u;
        auto const nyidx = i * 3u + 1u;
        auto const nzidx = i * 3u + 2u;

        auto const nx = static_cast<double>(geometry.normals[nxidx]);
        auto const ny = static_cast<double>(geometry.normals[nyidx]);
        auto const nz = static_cast<double>(geometry.normals[nzidx]);

        normals_.col(i) = Eigen::Vector3d{nx, ny, nz};
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

    auto const num_triangles = geometry.indices.size() / 3u;
    triangles_.resize(3u, num_triangles);
    for (std::size_t f = 0u; f < num_triangles; ++f)
    {
        auto const v1idx = f * 3u;
        auto const v2idx = f * 3u + 1u;
        auto const v3idx = f * 3u + 2u;

        auto const v1 = static_cast<std::uint32_t>(geometry.indices[v1idx]);
        auto const v2 = static_cast<std::uint32_t>(geometry.indices[v2idx]);
        auto const v3 = static_cast<std::uint32_t>(geometry.indices[v3idx]);

        triangles_.col(f) = triangle_type{v1, v2, v3};
    }
}

void shared_vertex_surface_mesh_t::set_color(Eigen::Vector3f const rgb)
{
    colors_.resizeLike(vertices_);
    colors_.colwise() = rgb;
}

void shared_vertex_surface_mesh_t::extract_normals()
{
    normals_.resizeLike(vertices_);
    normals_.setZero();

    index_type const num_boundary_faces = static_cast<index_type>(triangles_.cols());
    for (std::size_t f = 0u; f < num_boundary_faces; ++f)
    {
        triangle_type const triangle = triangles_.col(f);

        auto const v1 = triangle(0u);
        auto const v2 = triangle(1u);
        auto const v3 = triangle(2u);

        Eigen::Vector3d const p1 = vertices_.col(v1);
        Eigen::Vector3d const p2 = vertices_.col(v2);
        Eigen::Vector3d const p3 = vertices_.col(v3);

        Eigen::Vector3d const p21 = p2 - p1;
        Eigen::Vector3d const p31 = p3 - p1;

        Eigen::Vector3d const area_weighted_normal = 0.5 * p21.cross(p31);
        normals_.col(v1) += area_weighted_normal;
        normals_.col(v2) += area_weighted_normal;
        normals_.col(v3) += area_weighted_normal;
    }

    normals_.colwise().normalize();
}

std::vector<shared_vertex_surface_mesh_t::edge_type>
shared_vertex_surface_mesh_t::boundary_edges() const
{
    std::vector<edge_type> edges_with_duplicates{};

    std::size_t const num_triangles              = static_cast<std::size_t>(triangles_.cols());
    std::size_t constexpr num_edges_per_triangle = 3u;
    edges_with_duplicates.reserve(num_triangles * num_edges_per_triangle);

    std::map<edge_type, std::size_t> edge_occurrences{};

    for (std::size_t f = 0u; f < num_triangles; ++f)
    {
        triangle_type const triangle = triangles_.col(f);

        auto const v1 = triangle(0u);
        auto const v2 = triangle(1u);
        auto const v3 = triangle(2u);

        /**
         * Add index pairs in sorted order
         */
        auto const e1 = v1 < v2 ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
        auto const e2 = v2 < v3 ? std::make_pair(v2, v3) : std::make_pair(v3, v2);
        auto const e3 = v3 < v1 ? std::make_pair(v3, v1) : std::make_pair(v1, v3);

        ++edge_occurrences[e1];
        ++edge_occurrences[e2];
        ++edge_occurrences[e3];
    }

    std::vector<edge_type> boundary_edges{};
    boundary_edges.reserve(edge_occurrences.size()); // Over reserve memory
    for (auto const& [key, value] : edge_occurrences)
    {
        if (value > 1u)
            continue;

        boundary_edges.push_back(key);
    }

    return boundary_edges;
}

shared_vertex_surface_mesh_t::vertices_type const& shared_vertex_surface_mesh_t::vertices() const
{
    return vertices_;
}

shared_vertex_surface_mesh_t::vertices_type& shared_vertex_surface_mesh_t::vertices()
{
    return vertices_;
}

shared_vertex_surface_mesh_t::normals_type const& shared_vertex_surface_mesh_t::normals() const
{
    return normals_;
}

shared_vertex_surface_mesh_t::normals_type& shared_vertex_surface_mesh_t::normals()
{
    return normals_;
}

shared_vertex_surface_mesh_t::colors_type const& shared_vertex_surface_mesh_t::colors() const
{
    return colors_;
}

shared_vertex_surface_mesh_t::colors_type& shared_vertex_surface_mesh_t::colors()
{
    return colors_;
}

shared_vertex_surface_mesh_t::triangles_type const& shared_vertex_surface_mesh_t::triangles() const
{
    return triangles_;
}

shared_vertex_surface_mesh_t::triangles_type& shared_vertex_surface_mesh_t::triangles()
{
    return triangles_;
}

shared_vertex_surface_mesh_t::index_map_type const& shared_vertex_surface_mesh_t::index_map() const
{
    return index_map_;
}

shared_vertex_surface_mesh_t::index_map_type& shared_vertex_surface_mesh_t::index_map()
{
    return index_map_;
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

Eigen::Matrix3Xd rescale(
    Eigen::Matrix3Xd const& positions,
    Eigen::Vector3d const& boxmin,
    Eigen::Vector3d const& boxmax)
{
    Eigen::Vector3d const min = positions.rowwise().minCoeff();
    Eigen::Vector3d const max = positions.rowwise().maxCoeff();

    double const dx = max.x() - min.x();
    double const dy = max.y() - min.y();
    double const dz = max.z() - min.z();

    double constexpr eps           = 1e-8;
    bool const dx_division_by_zero = std::abs(dx) < eps;
    bool const dy_division_by_zero = std::abs(dy) < eps;
    bool const dz_division_by_zero = std::abs(dz) < eps;

    std::size_t const num_positions = static_cast<std::size_t>(positions.cols());

    Eigen::Matrix3Xd rescaled_positions(3u, num_positions);
    for (std::size_t i = 0u; i < num_positions; ++i)
    {
        Eigen::Vector3d const p = positions.col(i);

        auto const x = dx_division_by_zero ?
                           p.x() :
                           boxmin.x() + (boxmax.x() - boxmin.x()) * (p.x() - min.x()) / dx;
        auto const y = dy_division_by_zero ?
                           p.y() :
                           boxmin.y() + (boxmax.y() - boxmin.y()) * (p.y() - min.y()) / dy;
        auto const z = dz_division_by_zero ?
                           p.z() :
                           boxmin.z() + (boxmax.z() - boxmin.z()) * (p.z() - min.z()) / dz;

        rescaled_positions(0u, i) = x;
        rescaled_positions(1u, i) = y;
        rescaled_positions(2u, i) = z;
    }

    return rescaled_positions;
}

} // namespace common
} // namespace sbs