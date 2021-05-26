#include "geometry/get_simple_bar_model.h"

namespace sbs {
namespace geometry {

common::geometry_t get_simple_bar_model(std::size_t width, std::size_t height, std::size_t depth)
{
    common::geometry_t geometry;
    geometry.positions.resize(width * height * depth * 3u);

    for (auto i = 0; i < width; ++i)
    {
        for (auto j = 0; j < height; ++j)
        {
            for (auto k = 0; k < depth; ++k)
            {
                auto const row               = i * height * depth + j * depth + k;
                auto const idx               = row * 3u;
                geometry.positions[idx + 0u] = static_cast<float>(i);
                geometry.positions[idx + 1u] = static_cast<float>(j);
                geometry.positions[idx + 2u] = static_cast<float>(k);
            }
        }
    }

    auto const tet_count = (width - 1u) * (height - 1u) * (depth - 1u) * 5u;
    geometry.indices.resize(tet_count * 4u);
    for (std::size_t i = 0; i < width - 1u; i++)
    {
        for (std::size_t j = 0; j < height - 1u; j++)
        {
            for (std::size_t k = 0; k < depth - 1u; k++)
            {
                //     7*-----*6
                //     /|    /|
                //    / |   / |
                //  4*-----*5 |
                //   | 3*--|--*2
                //   | /   | /
                //   |/    |/
                //  0*-----*1

                // clang-format off
                int p0 = static_cast<int>(i        * height * depth + j        * depth + k       );
                int p1 = static_cast<int>((i + 1u) * height * depth + j        * depth + k       );
                int p2 = static_cast<int>((i + 1u) * height * depth + (j + 1u) * depth + k       );
                int p3 = static_cast<int>(i        * height * depth + (j + 1u) * depth + k       );
                int p4 = static_cast<int>(i        * height * depth + j        * depth + (k + 1u));
                int p5 = static_cast<int>((i + 1u) * height * depth + j        * depth + (k + 1u));
                int p6 = static_cast<int>((i + 1u) * height * depth + (j + 1u) * depth + (k + 1u));
                int p7 = static_cast<int>(i        * height * depth + (j + 1u) * depth + (k + 1u));
                // clang-format on

                auto const row = (i * (height - 1u) * (depth - 1u) + j * (depth - 1u) + k) * 5u;
                auto const tet_1_idx = row * 4u;
                auto const tet_2_idx = tet_1_idx + 1u * 4u;
                auto const tet_3_idx = tet_1_idx + 2u * 4u;
                auto const tet_4_idx = tet_1_idx + 3u * 4u;
                auto const tet_5_idx = tet_1_idx + 4u * 4u;

                if ((i + j + k) % 2 == 1)
                {
                    geometry.indices[tet_1_idx + 0u] = p1;
                    geometry.indices[tet_1_idx + 1u] = p0;
                    geometry.indices[tet_1_idx + 2u] = p5;
                    geometry.indices[tet_1_idx + 3u] = p2;

                    geometry.indices[tet_2_idx + 0u] = p5;
                    geometry.indices[tet_2_idx + 1u] = p2;
                    geometry.indices[tet_2_idx + 2u] = p7;
                    geometry.indices[tet_2_idx + 3u] = p6;

                    geometry.indices[tet_3_idx + 0u] = p7;
                    geometry.indices[tet_3_idx + 1u] = p0;
                    geometry.indices[tet_3_idx + 2u] = p5;
                    geometry.indices[tet_3_idx + 3u] = p4;

                    geometry.indices[tet_4_idx + 0u] = p2;
                    geometry.indices[tet_4_idx + 1u] = p0;
                    geometry.indices[tet_4_idx + 2u] = p7;
                    geometry.indices[tet_4_idx + 3u] = p3;

                    geometry.indices[tet_5_idx + 0u] = p5;
                    geometry.indices[tet_5_idx + 1u] = p0;
                    geometry.indices[tet_5_idx + 2u] = p7;
                    geometry.indices[tet_5_idx + 3u] = p2;
                }
                else
                {
                    geometry.indices[tet_1_idx + 0u] = p3;
                    geometry.indices[tet_1_idx + 1u] = p1;
                    geometry.indices[tet_1_idx + 2u] = p4;
                    geometry.indices[tet_1_idx + 3u] = p0;

                    geometry.indices[tet_2_idx + 0u] = p6;
                    geometry.indices[tet_2_idx + 1u] = p1;
                    geometry.indices[tet_2_idx + 2u] = p3;
                    geometry.indices[tet_2_idx + 3u] = p2;

                    geometry.indices[tet_3_idx + 0u] = p4;
                    geometry.indices[tet_3_idx + 1u] = p1;
                    geometry.indices[tet_3_idx + 2u] = p6;
                    geometry.indices[tet_3_idx + 3u] = p5;

                    geometry.indices[tet_4_idx + 0u] = p6;
                    geometry.indices[tet_4_idx + 1u] = p3;
                    geometry.indices[tet_4_idx + 2u] = p4;
                    geometry.indices[tet_4_idx + 3u] = p7;

                    geometry.indices[tet_5_idx + 0u] = p3;
                    geometry.indices[tet_5_idx + 1u] = p1;
                    geometry.indices[tet_5_idx + 2u] = p6;
                    geometry.indices[tet_5_idx + 3u] = p4;
                }
            }
        }
    }

    geometry.geometry_type = common::geometry_t::geometry_type_t::tetrahedron;

    return geometry;
}

} // namespace geometry
} // namespace sbs