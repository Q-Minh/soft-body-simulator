#ifndef SBS_COMMON_GEOMETRY_H
#define SBS_COMMON_GEOMETRY_H

#include "sbs/aliases.h"

#include <Eigen/Geometry>
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

    enum class geometry_type_t { triangle, tetrahedron, quad, hexahedron };
    geometry_type_t geometry_type;

    bool has_colors() const;
    bool has_positions() const;
    bool has_indices() const;
    bool has_normals() const;
    bool has_uvs() const;
    bool is_triangle_mesh() const;
    bool is_tetrahedral_mesh() const;

    void set_color(std::uint8_t r, std::uint8_t g, std::uint8_t b);
};

std::vector<Eigen::Vector3d> to_points(geometry_t const& geometry);

std::vector<index_type> to_indices(geometry_t const& geometry);

std::vector<Eigen::Vector3d> to_normals(geometry_t const& geometry);

std::vector<Eigen::Vector2f> to_uvs(geometry_t const& geometry);

std::vector<Eigen::Vector3f> to_colors(geometry_t const& geometry);

geometry_t transform(geometry_t const& geometry, Eigen::Affine3d const& transform);

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_GEOMETRY_H