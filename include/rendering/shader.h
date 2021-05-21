#ifndef SBS_RENDERING_SHADER_H
#define SBS_RENDERING_SHADER_H

#include <filesystem>
#include <fstream>
#include <glm/glm.hpp>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace sbs {
namespace rendering {

class shader_t
{
  public:
    shader_t() = default;

    shader_t(
        std::filesystem::path const& vertex_shader_path,
        std::filesystem::path const& fragment_shader_path);

    unsigned int id() const;
    bool should_use() const;
    void use() const;
    void set_bool_uniform(const std::string& name, bool value) const;
    void set_int_uniform(const std::string& name, int value) const;
    void set_float_uniform(const std::string& name, float value) const;
    void set_vec2_uniform(const std::string& name, const glm::vec2& value) const;
    void set_vec2_uniform(const std::string& name, float x, float y) const;
    void set_vec3_uniform(const std::string& name, const glm::vec3& value) const;
    void set_vec3_uniform(const std::string& name, float x, float y, float z) const;
    void set_vec4_uniform(const std::string& name, const glm::vec4& value) const;
    void set_vec4_uniform(const std::string& name, float x, float y, float z, float w);
    void set_mat2_uniform(const std::string& name, const glm::mat2& mat) const;
    void set_mat3_uniform(const std::string& name, const glm::mat3& mat) const;
    void set_mat4_uniform(const std::string& name, const glm::mat4& mat) const;
    void destroy() const;

  private:
    unsigned int id_;
    bool should_use_;
    std::vector<std::string> error_messages_;
};

} // namespace rendering
} // namespace sbs

#endif // SBS_RENDERING_SHADER_H