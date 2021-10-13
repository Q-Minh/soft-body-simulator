#include <sbs/geometry/get_simple_plane_model.h>

namespace sbs {
namespace geometry {

common::geometry_t get_simple_plane_model(
    std::array<scalar_type, 2u> const& min,
    std::array<scalar_type, 2u> const& max,
    scalar_type level)
{
    common::geometry_t geometry{};
    geometry.geometry_type = common::geometry_t::geometry_type_t::triangle;

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

    geometry.indices.push_back(0u);
    geometry.indices.push_back(1u);
    geometry.indices.push_back(3u);

    geometry.indices.push_back(1u);
    geometry.indices.push_back(2u);
    geometry.indices.push_back(3u);

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

    return geometry;
}

} // namespace geometry
} // namespace sbs