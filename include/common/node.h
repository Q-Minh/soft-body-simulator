#ifndef SBS_COMMON_NODE_H
#define SBS_COMMON_NODE_H

#include <vector>

namespace sbs {
namespace common {

/**
 * @brief
 */
struct node_t
{
    virtual void prepare_for_rendering() = 0;

    std::vector<float> positions;
    std::vector<unsigned int> indices;
    std::vector<float> normals;
    std::vector<float> uvs;
    std::vector<float> colors;

    enum class render_state_t { dirty, clean };
    render_state_t render_state;

    unsigned int VBO;
    unsigned int VAO;
    unsigned int EBO;
};

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_NODE_H