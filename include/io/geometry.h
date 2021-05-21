#ifndef SBS_IO_GEOMETRY_H
#define SBS_IO_GEOMETRY_H

#include <vector>

namespace sbs {
namespace io {

struct geometry_t
{
    std::vector<float> positions;
    std::vector<int> indices;
    std::vector<float> normals;
    std::vector<float> uvs;
    std::vector<std::uint8_t> colors;

    enum class geometry_type_t { triangle, tetrahedron };
    geometry_type_t geometry_type;
};

} // namespace io
} // namespace sbs

#endif // SBS_IO_GEOMETRY_H