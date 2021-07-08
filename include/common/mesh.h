#ifndef SBS_COMMON_MESH_H
#define SBS_COMMON_MESH_H

#include "node.h"

#include <tuple>

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

    virtual vertex_type vertex(std::size_t vi) const = 0;
    virtual triangle_type triangle(std::size_t f)    = 0;
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
    virtual triangle_type triangle(std::size_t f) override;
};

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_MESH_H