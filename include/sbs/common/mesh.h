#ifndef SBS_COMMON_MESH_H
#define SBS_COMMON_MESH_H

#include "node.h"

#include <Eigen/Core>
#include <array>

namespace sbs {
namespace common {

struct geometry_t;

class shared_vertex_surface_mesh_i : public renderable_node_t
{
  public:
    struct vertex_type
    {
        Eigen::Vector3d position;
        Eigen::Vector3d normal;
        Eigen::Vector3f color;
    };
    struct triangle_type
    {
        std::array<std::uint32_t, 3u> vertices;
    };

    virtual std::size_t triangle_count() const = 0;
    virtual std::size_t vertex_count() const   = 0;

    virtual vertex_type vertex(std::size_t vi) const    = 0;
    virtual triangle_type triangle(std::size_t f) const = 0;

    virtual void prepare_vertices_for_rendering() = 0;
    virtual void prepare_indices_for_rendering()  = 0;
};

class static_mesh : public shared_vertex_surface_mesh_i
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
 */
class dynamic_surface_mesh : public shared_vertex_surface_mesh_i
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

    vertex_type& mutable_vertex(std::size_t vi);
    triangle_type& mutable_triangle(std::size_t fi);

  private:
    std::vector<vertex_type> vertices_;
    std::vector<triangle_type> triangles_;
};

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_MESH_H