#include <sbs/geometry/get_simple_plane_model.h>

namespace sbs {
namespace geometry {

common::geometry_t get_simple_plane_model(
    std::array<scalar_type, 2u> const& min,
    std::array<scalar_type, 2u> const& max,
    scalar_type level,
    scalar_type thickness)
{
    common::geometry_t geometry{};
    geometry.geometry_type = common::geometry_t::geometry_type_t::triangle;

    /**
     *     y O
     *       |
     *       |
     *       O------O x
     *      /
     *     O
     *    z
     *
     *                            2                 1
     *                            o-----------------o
     *                        --- |             --- |
     *                    ---     |         ---     |
     *                 3 o-----------------o 0      |
     *                   |        |        |        |
     *                   |        |        |        |
     *                   |        |        |        |
     *                   |      6 o-----------------o 5
     *                   |    ---          |    ---
     *                   |---              |---
     *                 7 o-----------------o 4
     *
     */

    // top layer
    geometry.positions.push_back(static_cast<float>(max[0]));
    geometry.positions.push_back(static_cast<float>(level));
    geometry.positions.push_back(static_cast<float>(max[1]));

    geometry.positions.push_back(static_cast<float>(max[0]));
    geometry.positions.push_back(static_cast<float>(level));
    geometry.positions.push_back(static_cast<float>(min[1]));

    geometry.positions.push_back(static_cast<float>(min[0]));
    geometry.positions.push_back(static_cast<float>(level));
    geometry.positions.push_back(static_cast<float>(min[1]));

    geometry.positions.push_back(static_cast<float>(min[0]));
    geometry.positions.push_back(static_cast<float>(level));
    geometry.positions.push_back(static_cast<float>(max[1]));

    // bottom layer
    geometry.positions.push_back(static_cast<float>(max[0]));
    geometry.positions.push_back(static_cast<float>(level - thickness));
    geometry.positions.push_back(static_cast<float>(max[1]));

    geometry.positions.push_back(static_cast<float>(max[0]));
    geometry.positions.push_back(static_cast<float>(level - thickness));
    geometry.positions.push_back(static_cast<float>(min[1]));

    geometry.positions.push_back(static_cast<float>(min[0]));
    geometry.positions.push_back(static_cast<float>(level - thickness));
    geometry.positions.push_back(static_cast<float>(min[1]));

    geometry.positions.push_back(static_cast<float>(min[0]));
    geometry.positions.push_back(static_cast<float>(level - thickness));
    geometry.positions.push_back(static_cast<float>(max[1]));

    // triangles

    // top layer triangles
    geometry.indices.push_back(0u);
    geometry.indices.push_back(1u);
    geometry.indices.push_back(3u);

    geometry.indices.push_back(1u);
    geometry.indices.push_back(2u);
    geometry.indices.push_back(3u);

    // bottom layer triangles
    geometry.indices.push_back(4u);
    geometry.indices.push_back(7u);
    geometry.indices.push_back(5u);

    geometry.indices.push_back(5u);
    geometry.indices.push_back(7u);
    geometry.indices.push_back(6u);

    // right layer triangles
    geometry.indices.push_back(5u);
    geometry.indices.push_back(0u);
    geometry.indices.push_back(4u);

    geometry.indices.push_back(5u);
    geometry.indices.push_back(1u);
    geometry.indices.push_back(0u);

    // left layer triangles
    geometry.indices.push_back(7u);
    geometry.indices.push_back(2u);
    geometry.indices.push_back(6u);

    geometry.indices.push_back(7u);
    geometry.indices.push_back(3u);
    geometry.indices.push_back(2u);

    // back layer triangles
    geometry.indices.push_back(6u);
    geometry.indices.push_back(1u);
    geometry.indices.push_back(5u);

    geometry.indices.push_back(6u);
    geometry.indices.push_back(2u);
    geometry.indices.push_back(1u);

    // front layer triangles
    geometry.indices.push_back(4u);
    geometry.indices.push_back(3u);
    geometry.indices.push_back(7u);

    geometry.indices.push_back(4u);
    geometry.indices.push_back(0u);
    geometry.indices.push_back(3u);

    // top layer normals
    geometry.normals.push_back(0.f);
    geometry.normals.push_back(1.f);
    geometry.normals.push_back(0.f);

    geometry.normals.push_back(0.f);
    geometry.normals.push_back(1.f);
    geometry.normals.push_back(0.f);

    geometry.normals.push_back(0.f);
    geometry.normals.push_back(1.f);
    geometry.normals.push_back(0.f);

    geometry.normals.push_back(0.f);
    geometry.normals.push_back(1.f);
    geometry.normals.push_back(0.f);

    // bottom layer normals
    geometry.normals.push_back(0.f);
    geometry.normals.push_back(-1.f);
    geometry.normals.push_back(0.f);

    geometry.normals.push_back(0.f);
    geometry.normals.push_back(-1.f);
    geometry.normals.push_back(0.f);

    geometry.normals.push_back(0.f);
    geometry.normals.push_back(-1.f);
    geometry.normals.push_back(0.f);

    geometry.normals.push_back(0.f);
    geometry.normals.push_back(-1.f);
    geometry.normals.push_back(0.f);

    return geometry;
}

} // namespace geometry
} // namespace sbs