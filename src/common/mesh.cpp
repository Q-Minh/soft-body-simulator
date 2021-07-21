#include "sbs/common/mesh.h"

#include "sbs/common/geometry.h"
#include "sbs/common/primitive.h"

#include <array>
#include <numeric>

namespace sbs {
namespace common {

static_mesh::static_mesh(common::geometry_t const& geometry)
{
    bool const is_triangle_mesh =
        geometry.geometry_type == common::geometry_t::geometry_type_t::triangle;
    assert(is_triangle_mesh && "Should be surface triangle mesh");

    if (!is_triangle_mesh)
        return;

    assert(geometry.positions.size() % 3 == 0 && "must have (x,y,z) per position");
    assert(geometry.indices.size() % 3 == 0 && "must have (v1,v2,v3) per triangle");

    auto const vertex_count = geometry.indices.size();
    std::vector<float> vertex_buffer{};
    vertex_buffer.reserve(vertex_count * 9u);

    bool const has_colors = geometry.has_colors();

    assert(has_colors && "geometry for static mesh must have colors");

    if (has_colors)
    {
        assert(
            geometry.colors.size() == geometry.positions.size() &&
            "must have as many colors as positions");
    }

    for (std::size_t i = 0u; i < geometry.indices.size(); i += 3u)
    {
        auto const v1 = geometry.indices[i];
        auto const v2 = geometry.indices[i + 1u];
        auto const v3 = geometry.indices[i + 2u];

        auto const p1 = Eigen::Vector3d{
            geometry.positions[v1 * 3u],
            geometry.positions[v1 * 3u + 1u],
            geometry.positions[v1 * 3u + 2u]};
        auto const p2 = Eigen::Vector3d{
            geometry.positions[v2 * 3u],
            geometry.positions[v2 * 3u + 1u],
            geometry.positions[v2 * 3u + 2u]};
        auto const p3 = Eigen::Vector3d{
            geometry.positions[v3 * 3u],
            geometry.positions[v3 * 3u + 1u],
            geometry.positions[v3 * 3u + 2u]};

        auto const c1 = Eigen::Vector3f{
            static_cast<float>(geometry.colors[v1 * 3u]) / 255.f,
            static_cast<float>(geometry.colors[v1 * 3u + 1u]) / 255.f,
            static_cast<float>(geometry.colors[v1 * 3u + 2u]) / 255.f};
        auto const c2 = Eigen::Vector3f{
            static_cast<float>(geometry.colors[v2 * 3u]) / 255.f,
            static_cast<float>(geometry.colors[v2 * 3u + 1u]) / 255.f,
            static_cast<float>(geometry.colors[v2 * 3u + 2u]) / 255.f};
        auto const c3 = Eigen::Vector3f{
            static_cast<float>(geometry.colors[v3 * 3u]) / 255.f,
            static_cast<float>(geometry.colors[v3 * 3u + 1u]) / 255.f,
            static_cast<float>(geometry.colors[v3 * 3u + 2u]) / 255.f};

        common::triangle_t triangle_primitive{p1, p2, p3};

        auto const normal = triangle_primitive.normal();

        std::array<Eigen::Vector3f, 3u> const c{c1, c2, c3};
        std::array<Eigen::Vector3d, 3u> const p{p1, p2, p3};

        for (std::size_t j = 0u; j < p.size(); ++j)
        {
            vertex_buffer.push_back(static_cast<float>(p[j].x()));
            vertex_buffer.push_back(static_cast<float>(p[j].y()));
            vertex_buffer.push_back(static_cast<float>(p[j].z()));

            vertex_buffer.push_back(static_cast<float>(normal.x()));
            vertex_buffer.push_back(static_cast<float>(normal.y()));
            vertex_buffer.push_back(static_cast<float>(normal.z()));

            vertex_buffer.push_back(c[j].x());
            vertex_buffer.push_back(c[j].y());
            vertex_buffer.push_back(c[j].z());
        }
    }

    auto const index_count = geometry.indices.size();
    std::vector<std::uint32_t> index_buffer{};
    index_buffer.resize(index_count);
    std::iota(index_buffer.begin(), index_buffer.end(), 0u);

    transfer_vertices_for_rendering(std::move(vertex_buffer));
    transfer_indices_for_rendering(std::move(index_buffer));
}

void static_mesh::prepare_vertices_for_rendering()
{
    // no-op
}

void static_mesh::prepare_indices_for_rendering()
{
    // no-op
}

std::size_t static_mesh::triangle_count() const
{
    return get_cpu_index_buffer().size() / 3u;
}

std::size_t static_mesh::vertex_count() const
{
    return get_cpu_vertex_buffer().size() / 9u;
}

shared_vertex_surface_mesh_i::vertex_type static_mesh::vertex(std::size_t vi) const
{
    auto const& vertex_buffer = get_cpu_vertex_buffer();
    vertex_type v{};
    auto const idx = vi * 9u;
    v.x            = static_cast<double>(vertex_buffer[idx + 0u]);
    v.y            = static_cast<double>(vertex_buffer[idx + 1u]);
    v.z            = static_cast<double>(vertex_buffer[idx + 2u]);
    v.nx           = static_cast<double>(vertex_buffer[idx + 3u]);
    v.ny           = static_cast<double>(vertex_buffer[idx + 4u]);
    v.nz           = static_cast<double>(vertex_buffer[idx + 5u]);
    v.r            = vertex_buffer[idx + 6u];
    v.g            = vertex_buffer[idx + 7u];
    v.b            = vertex_buffer[idx + 8u];
    return v;
}

shared_vertex_surface_mesh_i::triangle_type static_mesh::triangle(std::size_t f) const
{
    auto const& index_buffer = get_cpu_index_buffer();
    triangle_type t{};
    auto const idx = f * 3u;
    t.v1           = index_buffer[idx + 0u];
    t.v2           = index_buffer[idx + 1u];
    t.v3           = index_buffer[idx + 2u];
    return t;
}

dynamic_surface_mesh::dynamic_surface_mesh(geometry_t const& geometry)
{
    auto const num_vertices = geometry.positions.size() / 3u;
    vertices_.reserve(num_vertices);

    assert(geometry.colors.size() == geometry.positions.size());

    auto const num_triangles = geometry.indices.size() / 3u;
    triangles_.reserve(num_triangles);

    for (std::size_t vi = 0u; vi < num_vertices; ++vi)
    {
        auto const idx = vi * 3u;

        vertex_type v{};
        v.x  = geometry.positions[idx + 0u];
        v.y  = geometry.positions[idx + 1u];
        v.z  = geometry.positions[idx + 2u];
        v.nx = geometry.has_normals() ? geometry.normals[idx + 0u] : 0.f;
        v.ny = geometry.has_normals() ? geometry.normals[idx + 1u] : 0.f;
        v.nz = geometry.has_normals() ? geometry.normals[idx + 2u] : 0.f;
        v.r  = static_cast<float>(geometry.colors[idx + 0u]) / 255.f;
        v.g  = static_cast<float>(geometry.colors[idx + 1u]) / 255.f;
        v.b  = static_cast<float>(geometry.colors[idx + 2u]) / 255.f;

        vertices_.push_back(v);
    }

    for (std::size_t fi = 0u; fi < num_triangles; ++fi)
    {
        auto const idx = fi * 3u;
        triangle_type f{};
        f.v1 = geometry.indices[idx + 0u];
        f.v2 = geometry.indices[idx + 1u];
        f.v3 = geometry.indices[idx + 2u];

        triangles_.push_back(f);
    }

    if (!geometry.has_normals())
    {
        for (auto const& triangle : triangles_)
        {
            auto& v1 = vertices_[triangle.v1];
            auto& v2 = vertices_[triangle.v2];
            auto& v3 = vertices_[triangle.v3];

            Eigen::Vector3d const p1{v1.x, v1.y, v1.z};
            Eigen::Vector3d const p2{v2.x, v2.y, v2.z};
            Eigen::Vector3d const p3{v3.x, v3.y, v3.z};

            common::triangle_t const triangle_primitive{p1, p2, p3};

            auto const n = triangle_primitive.normal();
            auto const A = triangle_primitive.area();

            v1.nx += A * n.x();
            v1.ny += A * n.y();
            v1.nz += A * n.z();
            v2.nx += A * n.x();
            v2.ny += A * n.y();
            v2.nz += A * n.z();
            v3.nx += A * n.x();
            v3.ny += A * n.y();
            v3.nz += A * n.z();
        }
        std::for_each(vertices_.begin(), vertices_.end(), [](vertex_type& v) {
            Eigen::Vector3d const n = Eigen::Vector3d{v.nx, v.ny, v.nz}.normalized();
            v.nx                    = n.x();
            v.ny                    = n.y();
            v.nz                    = n.z();
        });
    }
}

void dynamic_surface_mesh::prepare_vertices_for_rendering()
{
    auto const vertex_count = vertices_.size();
    std::vector<float> vertex_buffer{};
    vertex_buffer.reserve(vertex_count * 9u);

    for (auto const& v : vertices_)
    {
        vertex_buffer.push_back(static_cast<float>(v.x));
        vertex_buffer.push_back(static_cast<float>(v.y));
        vertex_buffer.push_back(static_cast<float>(v.z));
        vertex_buffer.push_back(static_cast<float>(v.nx));
        vertex_buffer.push_back(static_cast<float>(v.ny));
        vertex_buffer.push_back(static_cast<float>(v.nz));
        vertex_buffer.push_back(v.r);
        vertex_buffer.push_back(v.g);
        vertex_buffer.push_back(v.b);
    }

    transfer_vertices_for_rendering(std::move(vertex_buffer));
}

void dynamic_surface_mesh::prepare_indices_for_rendering()
{
    auto const index_count = triangles_.size() * 3u;
    std::vector<std::uint32_t> index_buffer{};
    index_buffer.reserve(index_count);

    for (auto const& triangle : triangles_)
    {
        index_buffer.push_back(triangle.v1);
        index_buffer.push_back(triangle.v2);
        index_buffer.push_back(triangle.v3);
    }

    transfer_indices_for_rendering(std::move(index_buffer));
}

std::size_t dynamic_surface_mesh::triangle_count() const
{
    return triangles_.size();
}

std::size_t dynamic_surface_mesh::vertex_count() const
{
    return vertices_.size();
}

dynamic_surface_mesh::vertex_type dynamic_surface_mesh::vertex(std::size_t vi) const
{
    return vertices_[vi];
}

dynamic_surface_mesh::triangle_type dynamic_surface_mesh::triangle(std::size_t fi) const
{
    return triangles_[fi];
}

dynamic_surface_mesh::vertex_type& dynamic_surface_mesh::mutable_vertex(std::size_t vi)
{
    return vertices_[vi];
}

dynamic_surface_mesh::triangle_type& dynamic_surface_mesh::mutable_triangle(std::size_t fi)
{
    return triangles_[fi];
}

} // namespace common
} // namespace sbs
