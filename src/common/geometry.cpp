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

} // namespace common
} // namespace sbs
