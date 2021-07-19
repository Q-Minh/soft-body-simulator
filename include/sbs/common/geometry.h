#ifndef SBS_COMMON_GEOMETRY_H
#define SBS_COMMON_GEOMETRY_H

#include <vector>

namespace sbs {
namespace common {

struct geometry_t
{
    std::vector<float> positions;
    std::vector<int> indices;
    std::vector<float> normals;
    std::vector<float> uvs;
    std::vector<std::uint8_t> colors;

    enum class geometry_type_t { triangle, tetrahedron };
    geometry_type_t geometry_type;

    bool has_colors() const;
    bool has_positions() const;
    bool has_indices() const;
    bool has_normals() const;
    bool has_uvs() const;
    bool is_triangle_mesh() const;
    bool is_tetrahedral_mesh() const;
};

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_GEOMETRY_H