#ifndef SBS_GEOMETRY_GET_SIMPLE_PLANE_MODEL_H
#define SBS_GEOMETRY_GET_SIMPLE_PLANE_MODEL_H

#include <array>
#include <sbs/aliases.h>
#include <sbs/common/geometry.h>

namespace sbs {
namespace geometry {

common::geometry_t get_simple_plane_model(
    std::array<scalar_type, 2u> const& min,
    std::array<scalar_type, 2u> const& max,
    scalar_type level, 
    scalar_type thickness);

} // namespace geometry
} // namespace sbs

#endif // SBS_GEOMETRY_GET_SIMPLE_PLANE_MODEL_H