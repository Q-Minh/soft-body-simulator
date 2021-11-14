#ifndef SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_MLS_SURFACE_H
#define SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_MLS_SURFACE_H

#include <Eigen/Core>
#include <array>
#include <optional>
#include <sbs/common/mesh.h>
#include <vector>

namespace sbs {

namespace topology {
class triangle_t;
} // namespace topology

namespace physics {
namespace mechanics {

class hybrid_mesh_meshless_mls_body_t;
class hybrid_mesh_meshless_mls_node_t;

class hybrid_mesh_meshless_mls_surface_vertex_t
{
  public:
    hybrid_mesh_meshless_mls_surface_vertex_t(Eigen::Vector3d const& x0)
        : ti_(std::numeric_limits<index_type>::max()), Xi_(x0), xi_(x0)
    {
    }

    void initialize(
        std::vector<index_type> const& meshless_neighbours,
        hybrid_mesh_meshless_mls_body_t const& mechanical_model,
        index_type ti);

    index_type ti() const { return ti_; }
    index_type& ti() { return ti_; }

    bool is_in_tetrahedron() const { return ti_ != std::numeric_limits<index_type>::max(); }

    Eigen::Vector3d const& xi() const;
    Eigen::Vector3d& xi();
    Eigen::Vector3d const& Xi() const;
    Eigen::Vector3d& Xi();

    std::vector<scalar_type> const& phi_js() const;
    std::array<std::optional<scalar_type>, 4u> const& mesh_phi_js() const;

    std::vector<index_type> const& neighbours() const;

  private:
    index_type ti_;
    Eigen::Vector3d Xi_;
    Eigen::Vector3d xi_;
    std::vector<scalar_type> phi_js_;
    std::array<std::optional<scalar_type>, 4u> mesh_phi_js_;
    std::vector<index_type> neighbours_;
};

class hybrid_mesh_meshless_sph_surface_t : public common::shared_vertex_surface_mesh_i
{
  public:
    using vertex_type   = common::shared_vertex_surface_mesh_i::vertex_type;
    using triangle_type = common::shared_vertex_surface_mesh_i::triangle_type;

    hybrid_mesh_meshless_sph_surface_t() = default;
    hybrid_mesh_meshless_sph_surface_t(
        hybrid_mesh_meshless_mls_body_t* mechanical_model,
        std::vector<Eigen::Vector3d> const& vertices,
        std::vector<topology::triangle_t> const& triangles);

    hybrid_mesh_meshless_sph_surface_t(hybrid_mesh_meshless_sph_surface_t const& other) = default;
    hybrid_mesh_meshless_sph_surface_t(hybrid_mesh_meshless_sph_surface_t&& other)      = default;
    hybrid_mesh_meshless_sph_surface_t&
    operator=(hybrid_mesh_meshless_sph_surface_t const& other) = default;
    hybrid_mesh_meshless_sph_surface_t&
    operator=(hybrid_mesh_meshless_sph_surface_t&& other) = default;

    virtual std::size_t vertex_count() const override;
    virtual std::size_t triangle_count() const override;

    virtual vertex_type vertex(std::size_t vi) const override;
    virtual triangle_type triangle(std::size_t f) const override;

    virtual void prepare_vertices_for_rendering() override;
    virtual void prepare_indices_for_rendering() override;

    void initialize_interpolation_scheme(scalar_type const h);

    vertex_type& world_space_vertex(std::size_t vi);
    Eigen::Vector3d& material_space_position(std::size_t vi);

    std::vector<hybrid_mesh_meshless_mls_surface_vertex_t> const& embedded_vertices() const;

    hybrid_mesh_meshless_mls_body_t const* mechanical_model() const;
    hybrid_mesh_meshless_mls_body_t* mechanical_model();

    void compute_positions();
    void compute_normals();

  private:
    void prepare_vertices_for_surface_rendering();

    std::vector<vertex_type>
        render_vertices_; ///< The vertices_ member of tetrahedral_mesh_boundary_t will be
                          ///< considered as vertex positions in material space, while the
                          ///< render_vertices_ are the world space vertex positions
    std::vector<triangle_type> triangles_;
    hybrid_mesh_meshless_mls_body_t* mechanical_model_;

    std::vector<hybrid_mesh_meshless_mls_surface_vertex_t> vertices_;
};

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_MLS_SURFACE_H