#include "geometry/get_simple_cloth_model.h"

namespace sbs {
namespace geometry {

common::geometry_t get_simple_cloth_model(int rows, int cols)
{
    common::geometry_t geometry;

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            auto const xoffset = static_cast<float>(i);
            auto const yoffset = static_cast<float>(j);

            geometry.positions.push_back(xoffset);
            geometry.positions.push_back(yoffset);
            geometry.positions.push_back(0.);

            if (i == rows - 1u || j == cols - 1u)
                continue;

            auto const lower_left_corner  = i * cols + j;
            auto const upper_left_corner  = i * cols + (j + 1);
            auto const lower_right_corner = (i + 1) * cols + j;
            auto const upper_right_corner = (i + 1) * cols + (j + 1);

            geometry.indices.push_back(lower_left_corner);
            geometry.indices.push_back(upper_right_corner);
            geometry.indices.push_back(upper_left_corner);

            geometry.indices.push_back(lower_left_corner);
            geometry.indices.push_back(lower_right_corner);
            geometry.indices.push_back(upper_right_corner);
        }
    }

    geometry.geometry_type = common::geometry_t::geometry_type_t::triangle;

    return geometry;
}

} // namespace geometry
} // namespace sbs