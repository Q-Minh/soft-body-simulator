#ifndef SBS_COMMON_MESH_H
#define SBS_COMMON_MESH_H

#include "node.h"

namespace sbs {
namespace common {

struct geometry_t;

class shared_vertex_surface_mesh_i
{
  public:
    struct vertex_type
    {
        double x, y, z;
        double nx, ny, nz;
        float r, g, b;
    };
    struct triangle_type
    {
        std::uint32_t v1, v2, v3;
    };

    virtual std::size_t triangle_count() const = 0;
    virtual std::size_t vertex_count() const   = 0;

    virtual vertex_type vertex(std::size_t vi) const    = 0;
    virtual triangle_type triangle(std::size_t f) const = 0;
};

class static_mesh : public renderable_node_t, public shared_vertex_surface_mesh_i
{
  public:
    static_mesh() = default;
    static_mesh(geometry_t const& geometry);

    virtual void prepare_vertices_for_rendering() override;
    virtual void prepare_indices_for_rendering() override;

    virtual std::size_t triangle_count() const override;
    virtual std::size_t vertex_count() const override;

    virtual vertex_type vertex(std::size_t vi) const override;
    virtual triangle_type triangle(std::size_t fi) const override;
};

/**
 * @brief Shared vertex surface mesh that can have its vertices and indices modified and rendered.
 * TODO: Still not implemented
 */
class dynamic_surface_mesh : public renderable_node_t, public shared_vertex_surface_mesh_i
{
  public:
    using vertex_type   = shared_vertex_surface_mesh_i::vertex_type;
    using triangle_type = shared_vertex_surface_mesh_i::triangle_type;

    dynamic_surface_mesh() = default;
    dynamic_surface_mesh(geometry_t const& geometry);

    virtual void prepare_vertices_for_rendering() override;
    virtual void prepare_indices_for_rendering() override;

    virtual std::size_t triangle_count() const override;
    virtual std::size_t vertex_count() const override;

    virtual vertex_type vertex(std::size_t vi) const override;
    virtual triangle_type triangle(std::size_t fi) const override;

  private:
    std::vector<vertex_type> vertices_;
    std::vector<triangle_type> triangles_;
};

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_MESH_H