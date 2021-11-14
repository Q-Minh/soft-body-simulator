#include "sbs/common/geometry.h"

namespace sbs {
namespace common {

bool geometry_t::has_colors() const
{
    return !colors.empty();
}

bool geometry_t::has_positions() const
{
    return !positions.empty();
}

bool geometry_t::has_indices() const
{
    return !indices.empty();
}

bool geometry_t::has_normals() const
{
    return !normals.empty();
}

bool geometry_t::has_uvs() const
{
    return !uvs.empty();
}

bool geometry_t::is_triangle_mesh() const
{
    return geometry_type == geometry_type_t::triangle;
}

bool geometry_t::is_tetrahedral_mesh() const
{
    return geometry_type == geometry_type_t::tetrahedron;
}

void geometry_t::set_color(std::uint8_t r, std::uint8_t g, std::uint8_t b)
{
    colors.clear();
    for (std::size_t i = 0u; i < positions.size(); i += 3u)
    {
        colors.push_back(r);
        colors.push_back(g);
        colors.push_back(b);
    }
}

std::vector<Eigen::Vector3d> to_points(geometry_t const& geometry)
{
    std::vector<Eigen::Vector3d> points{};
    points.reserve(geometry.positions.size() / 3u);
    for (auto i = 0u; i < geometry.positions.size(); i += 3u)
    {
        Eigen::Vector3d const p{
            geometry.positions[i],
            geometry.positions[i + 1],
            geometry.positions[i + 2]};
        points.push_back(p);
    }
    return points;
}

std::vector<index_type> to_indices(geometry_t const& geometry)
{
    std::vector<index_type> indices{};
    indices.reserve(geometry.indices.size());
    for (auto idx : geometry.indices)
    {
        indices.push_back(static_cast<index_type>(idx));
    }
    return indices;
}

std::vector<Eigen::Vector3d> to_normals(geometry_t const& geometry)
{
    std::vector<Eigen::Vector3d> normals{};
    normals.reserve(geometry.normals.size() / 3u);
    for (auto i = 0u; i < geometry.normals.size(); i += 3u)
    {
        Eigen::Vector3d const n{
            geometry.normals[i],
            geometry.normals[i + 1],
            geometry.normals[i + 2]};
        normals.push_back(n);
    }
    return normals;
}

std::vector<Eigen::Vector2f> to_uvs(geometry_t const& geometry)
{
    std::vector<Eigen::Vector2f> uvs{};
    uvs.reserve(geometry.uvs.size() / 2u);
    for (auto i = 0u; i < geometry.uvs.size(); i += 2u)
    {
        Eigen::Vector2f const uv{geometry.uvs[i], geometry.uvs[i + 1]};
        uvs.push_back(uv);
    }
    return uvs;
}

std::vector<Eigen::Vector3f> to_colors(geometry_t const& geometry)
{
    std::vector<Eigen::Vector3f> colors{};
    colors.reserve(geometry.colors.size() / 3u);
    for (auto i = 0u; i < geometry.colors.size(); i += 3u)
    {
        Eigen::Vector3f const c{
            static_cast<float>(geometry.colors[i]) / 255.f,
            static_cast<float>(geometry.colors[i + 1]) / 255.f,
            static_cast<float>(geometry.colors[i + 2]) / 255.f};
        colors.push_back(c);
    }
    return colors;
}

geometry_t transform(geometry_t const& geometry, Eigen::Affine3d const& transform)
{
    geometry_t other = geometry;
    for (auto i = 0u; i < geometry.positions.size(); i += 3u)
    {
        Eigen::Vector3d const p{
            geometry.positions[i],
            geometry.positions[i + 1],
            geometry.positions[i + 2]};

        Eigen::Vector3d const p_prime = transform * p;
        other.positions[i]            = p_prime.x();
        other.positions[i + 1]        = p_prime.y();
        other.positions[i + 2]        = p_prime.z();
    }
    return other;
}

} // namespace common
} // namespace sbs
