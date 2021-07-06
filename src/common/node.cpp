#include "common/node.h"

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

std::vector<float> const& renderable_node_t::get_render_vertices() const
{
    return vertices_;
}

std::vector<std::uint32_t> const& renderable_node_t::get_render_indices() const
{
    return indices_;
}

renderable_node_t::position_t renderable_node_t::get_position(unsigned int v) const
{
    unsigned int const offset = 9u * v;
    position_t p{};
    p.x = vertices_[offset];
    p.y = vertices_[offset + 1u];
    p.z = vertices_[offset + 2u];
    return p;
}

renderable_node_t::normal_t renderable_node_t::get_normal(unsigned int v) const
{
    unsigned int const offset = 9u * v;
    normal_t n{};
    n.nx = vertices_[offset + 3u];
    n.ny = vertices_[offset + 4u];
    n.nz = vertices_[offset + 5u];
    return n;
}

renderable_node_t::triangle_t renderable_node_t::get_triangle(unsigned int t) const
{
    unsigned int const offset = 3u * t;
    triangle_t f{};
    f.v1 = indices_[offset + 0u];
    f.v2 = indices_[offset + 1u];
    f.v3 = indices_[offset + 2u];
    return f;
}

} // namespace common
} // namespace sbs