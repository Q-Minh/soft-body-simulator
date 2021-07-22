#include "sbs/rendering/shader.h"

#include <glad/glad.h>

namespace sbs {
namespace rendering {

const char* shader_t::vertex_shader_position_attribute_name = "VertexShaderVertexPosition";
const char* shader_t::vertex_shader_normal_attribute_name   = "VertexShaderVertexNormal";
const char* shader_t::vertex_shader_color_attribute_name    = "VertexShaderVertexColor";

shader_t::shader_t(
    std::filesystem::path const& vertex_shader_path,
    std::filesystem::path const& fragment_shader_path)
    : id_(0u), should_use_{false}
{
    bool const is_vertex_shader_path_valid =
        std::filesystem::exists(vertex_shader_path) && vertex_shader_path.has_filename() &&
        vertex_shader_path.has_extension() && vertex_shader_path.extension().string() == ".vs";

    bool const is_fragment_shader_path_valid =
        std::filesystem::exists(fragment_shader_path) && fragment_shader_path.has_filename() &&
        fragment_shader_path.has_extension() && fragment_shader_path.extension().string() == ".fs";

    if (!is_vertex_shader_path_valid || !is_fragment_shader_path_valid)
    {
        error_messages_.push_back("Invalid shader paths");
        return;
    }

    std::ifstream vertex_ifs{vertex_shader_path.string()};
    std::ifstream fragment_ifs{fragment_shader_path.string()};

    std::string const vertex_shader_source_code{
        std::istreambuf_iterator<char>(vertex_ifs),
        std::istreambuf_iterator<char>{}};
    std::string const fragment_shader_source_code{
        std::istreambuf_iterator<char>(fragment_ifs),
        std::istreambuf_iterator<char>{}};

    if (!vertex_ifs.is_open())
    {
        error_messages_.push_back("Cannot open vertex shader");
        return;
    }

    if (!fragment_ifs.is_open())
    {
        error_messages_.push_back("Cannot open fragment shader");
        return;
    }

    /**
     * Compile vertex shader
     */
    unsigned int vertex_shader_id              = glCreateShader(GL_VERTEX_SHADER);
    const char* vertex_shader_source_code_cstr = vertex_shader_source_code.c_str();
    glShaderSource(vertex_shader_id, 1, &vertex_shader_source_code_cstr, NULL);
    glCompileShader(vertex_shader_id);

    GLint vertex_shader_compilation_success;
    char vertex_shader_compilation_info[512];
    glGetShaderiv(vertex_shader_id, GL_COMPILE_STATUS, &vertex_shader_compilation_success);
    if (!vertex_shader_compilation_success)
    {
        glGetShaderInfoLog(vertex_shader_id, 512, NULL, vertex_shader_compilation_info);
        std::ostringstream oss{};
        oss << "Vertex shader compilation error:\n" << vertex_shader_compilation_info << "\n";
        error_messages_.push_back(oss.str());
    }

    /**
     * Compile fragment shader
     */
    unsigned int fragment_shader_id              = glCreateShader(GL_FRAGMENT_SHADER);
    const char* fragment_shader_source_code_cstr = fragment_shader_source_code.c_str();
    glShaderSource(fragment_shader_id, 1, &fragment_shader_source_code_cstr, NULL);
    glCompileShader(fragment_shader_id);

    GLint fragment_shader_compilation_success;
    char fragment_shader_compilation_info[512];
    glGetShaderiv(fragment_shader_id, GL_COMPILE_STATUS, &fragment_shader_compilation_success);
    if (!fragment_shader_compilation_success)
    {
        glGetShaderInfoLog(fragment_shader_id, 512, nullptr, fragment_shader_compilation_info);
        std::ostringstream oss{};
        oss << "Fragment shader compilation error:\n" << fragment_shader_compilation_info << "\n";
        error_messages_.push_back(oss.str());
    }

    unsigned int shader_program_id = glCreateProgram();
    glAttachShader(shader_program_id, vertex_shader_id);
    glAttachShader(shader_program_id, fragment_shader_id);
    glLinkProgram(shader_program_id);

    GLint shader_program_link_success;
    char shader_program_link_info[512];
    glGetProgramiv(shader_program_id, GL_LINK_STATUS, &shader_program_link_success);
    if (!shader_program_link_success)
    {
        glGetProgramInfoLog(shader_program_id, 512, nullptr, shader_program_link_info);
        std::ostringstream oss{};
        oss << "Shader program linking error:\n" << shader_program_link_info << "\n";
        error_messages_.push_back(oss.str());
    }

    should_use_ = static_cast<bool>(vertex_shader_compilation_success) &&
                  static_cast<bool>(fragment_shader_compilation_success) &&
                  static_cast<bool>(shader_program_link_success);

    if (!should_use_)
    {
        return;
    }

    glDeleteShader(vertex_shader_id);
    glDeleteShader(fragment_shader_id);

    id_ = shader_program_id;
}

unsigned int shader_t::id() const
{
    return id_;
}

bool shader_t::should_use() const
{
    return should_use_;
}

void shader_t::use() const
{
    glUseProgram(id_);
}

void shader_t::set_bool_uniform(const std::string& name, bool value) const
{
    glUniform1i(glGetUniformLocation(id_, name.c_str()), (int)value);
}

void shader_t::set_int_uniform(const std::string& name, int value) const
{
    glUniform1i(glGetUniformLocation(id_, name.c_str()), value);
}

void shader_t::set_float_uniform(const std::string& name, float value) const
{
    glUniform1f(glGetUniformLocation(id_, name.c_str()), value);
}

void shader_t::set_vec2_uniform(const std::string& name, const glm::vec2& value) const
{
    glUniform2fv(glGetUniformLocation(id_, name.c_str()), 1, &value[0]);
}

void shader_t::set_vec2_uniform(const std::string& name, float x, float y) const
{
    glUniform2f(glGetUniformLocation(id_, name.c_str()), x, y);
}

void shader_t::set_vec3_uniform(const std::string& name, const glm::vec3& value) const
{
    glUniform3fv(glGetUniformLocation(id_, name.c_str()), 1, &value[0]);
}

void shader_t::set_vec3_uniform(const std::string& name, float x, float y, float z) const
{
    glUniform3f(glGetUniformLocation(id_, name.c_str()), x, y, z);
}

void shader_t::set_vec4_uniform(const std::string& name, const glm::vec4& value) const
{
    glUniform4fv(glGetUniformLocation(id_, name.c_str()), 1, &value[0]);
}

void shader_t::set_vec4_uniform(const std::string& name, float x, float y, float z, float w)
{
    glUniform4f(glGetUniformLocation(id_, name.c_str()), x, y, z, w);
}

void shader_t::set_mat2_uniform(const std::string& name, const glm::mat2& mat) const
{
    glUniformMatrix2fv(glGetUniformLocation(id_, name.c_str()), 1, GL_FALSE, &mat[0][0]);
}

void shader_t::set_mat3_uniform(const std::string& name, const glm::mat3& mat) const
{
    glUniformMatrix3fv(glGetUniformLocation(id_, name.c_str()), 1, GL_FALSE, &mat[0][0]);
}

void shader_t::set_mat4_uniform(const std::string& name, const glm::mat4& mat) const
{
    glUniformMatrix4fv(glGetUniformLocation(id_, name.c_str()), 1, GL_FALSE, &mat[0][0]);
}

void shader_t::destroy() const
{
    glDeleteProgram(id_);
}

} // namespace rendering
} // namespace sbs