#include "geometry.h"

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

} // namespace common
} // namespace sbs
