#include "sbs/common/node.h"

namespace sbs {
namespace common {

void renderable_node_t::set_id(std::string const& id)
{
    id_ = id;
}

std::string const& renderable_node_t::id() const
{
    return id_;
}

void renderable_node_t::set_vao(unsigned int vao)
{
    VAO_ = vao;
}

void renderable_node_t::set_vbo(unsigned int vbo)
{
    VBO_ = vbo;
}

void renderable_node_t::set_ebo(unsigned int ebo)
{
    EBO_ = ebo;
}

unsigned int const& renderable_node_t::VAO() const
{
    return VAO_;
}

unsigned int const& renderable_node_t::VBO() const
{
    return VBO_;
}

unsigned int const& renderable_node_t::EBO() const
{
    return EBO_;
}

unsigned int& renderable_node_t::VAO()
{
    return VAO_;
}

unsigned int& renderable_node_t::VBO()
{
    return VBO_;
}

unsigned int& renderable_node_t::EBO()
{
    return EBO_;
}

void renderable_node_t::mark_vertices_dirty()
{
    render_state_.should_transfer_vertices = true;
}

void renderable_node_t::mark_indices_dirty()
{
    render_state_.should_transfer_indices = true;
}

void renderable_node_t::mark_should_render_wireframe()
{
    render_state_.should_render_wireframe = true;
}

void renderable_node_t::mark_vertices_clean()
{
    render_state_.should_transfer_vertices = false;
}

void renderable_node_t::mark_indices_clean()
{
    render_state_.should_transfer_indices = false;
}

void renderable_node_t::mark_should_render_triangles()
{
    render_state_.should_render_wireframe = false;
}

bool renderable_node_t::should_transfer_vertices() const
{
    return render_state_.should_transfer_vertices;
}

bool renderable_node_t::should_transfer_indices() const
{
    return render_state_.should_transfer_indices;
}

bool renderable_node_t::should_render_wireframe() const
{
    return render_state_.should_render_wireframe;
}

bool renderable_node_t::is_environment_body() const
{
    return body_type_ == body_type_t::environment;
}

bool renderable_node_t::is_physically_simulated_body() const
{
    return body_type_ == body_type_t::physical;
}

void renderable_node_t::set_as_environment_body()
{
    body_type_ = body_type_t::environment;
}

void renderable_node_t::set_as_physically_simulated_body()
{
    body_type_ = body_type_t::physical;
}

void renderable_node_t::set_as_collideable_body()
{
    is_collideable_ = true;
}

void renderable_node_t::set_as_non_collideable_body()
{
    is_collideable_ = false;
}

bool renderable_node_t::is_collideable_body() const
{
    return is_collideable_;
}

std::vector<float> const& renderable_node_t::get_cpu_vertex_buffer() const
{
    return cpu_vertex_buffer_;
}

std::vector<std::uint32_t> const& renderable_node_t::get_cpu_index_buffer() const
{
    return cpu_index_buffer_;
}

void renderable_node_t::transfer_vertices_for_rendering(std::vector<float>&& vertices)
{
    cpu_vertex_buffer_ = std::move(vertices);
}

void renderable_node_t::transfer_indices_for_rendering(std::vector<std::uint32_t>&& indices)
{
    cpu_index_buffer_ = std::move(indices);
}

} // namespace common
} // namespace sbs