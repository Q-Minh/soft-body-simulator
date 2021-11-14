#ifndef SBS_PHYSICS_MECHANICS_MESHLESS_SURFACE_H
#define SBS_PHYSICS_MECHANICS_MESHLESS_SURFACE_H

#include <Eigen/Core>
#include <sbs/common/mesh.h>
#include <sbs/topology/tetrahedron_set.h>

namespace sbs {
namespace physics {
namespace mechanics {

class meshless_sph_body_t;

class meshless_sph_surface_vertex_t
{
  public:
    meshless_sph_surface_vertex_t(meshless_sph_surface_vertex_t const& other) = default;
    meshless_sph_surface_vertex_t(meshless_sph_surface_vertex_t&& other)      = default;

    meshless_sph_surface_vertex_t& operator=(meshless_sph_surface_vertex_t const& other) = default;
    meshless_sph_surface_vertex_t& operator=(meshless_sph_surface_vertex_t&& other) = default;

    meshless_sph_surface_vertex_t(Eigen::Vector3d const& x0)
        : x0_(x0), x_(), Xkjs_(), Wkjs_(), Vjs_(), sk_(), neighbours_()
    {
    }

    meshless_sph_surface_vertex_t(
        Eigen::Vector3d const& x0,
        Eigen::Vector3d const& x,
        std::vector<Eigen::Vector3d> const& Xkjs,
        std::vector<scalar_type> const& Wkjs,
        std::vector<scalar_type> const& Vjs,
        scalar_type const sk,
        std::vector<index_type> const& neighbours)
        : x0_(x0), x_(x), Xkjs_(Xkjs), Wkjs_(Wkjs), Vjs_(Vjs), sk_(sk), neighbours_(neighbours)
    {
    }

    Eigen::Vector3d const& x0() const { return x0_; }
    Eigen::Vector3d& x0() { return x0_; }
    Eigen::Vector3d const& x() const { return x_; }
    Eigen::Vector3d& x() { return x_; }
    std::vector<Eigen::Vector3d> const& Xkjs() const { return Xkjs_; }
    std::vector<Eigen::Vector3d>& Xkjs() { return Xkjs_; }
    std::vector<scalar_type> const& Wkjs() const { return Wkjs_; }
    std::vector<scalar_type>& Wkjs() { return Wkjs_; }
    std::vector<scalar_type> const& Vjs() const { return Vjs_; }
    std::vector<scalar_type>& Vjs() { return Vjs_; }
    scalar_type sk() const { return sk_; }
    scalar_type& sk() { return sk_; }
    std::vector<index_type> const& neighbours() const { return neighbours_; }
    std::vector<index_type>& neighbours() { return neighbours_; }

  private:
    Eigen::Vector3d x0_;
    Eigen::Vector3d x_;
    std::vector<Eigen::Vector3d> Xkjs_;
    std::vector<scalar_type> Wkjs_;
    std::vector<scalar_type> Vjs_;
    scalar_type sk_;
    std::vector<index_type> neighbours_;
};

class meshless_sph_surface_t : public common::shared_vertex_surface_mesh_i
{
  public:
    using vertex_type   = common::shared_vertex_surface_mesh_i::vertex_type;
    using triangle_type = common::shared_vertex_surface_mesh_i::triangle_type;

    meshless_sph_surface_t() = default;
    meshless_sph_surface_t(
        meshless_sph_body_t* mechanical_model,
        std::vector<Eigen::Vector3d> const& vertices,
        std::vector<topology::triangle_t> const& triangles);

    meshless_sph_surface_t(meshless_sph_surface_t const& other) = default;
    meshless_sph_surface_t(meshless_sph_surface_t&& other)      = default;
    meshless_sph_surface_t& operator=(meshless_sph_surface_t const& other) = default;
    meshless_sph_surface_t& operator=(meshless_sph_surface_t&& other) = default;

    void initialize_interpolation_scheme(scalar_type const h);

    virtual std::size_t triangle_count() const override;
    virtual std::size_t vertex_count() const override;

    virtual vertex_type vertex(std::size_t vi) const override;
    virtual triangle_type triangle(std::size_t f) const override;

    virtual void prepare_vertices_for_rendering() override;
    virtual void prepare_indices_for_rendering() override;

    vertex_type& world_space_vertex(std::size_t vi);
    Eigen::Vector3d& material_space_position(std::size_t vi);
    void compute_positions();
    void compute_normals();

    std::vector<meshless_sph_surface_vertex_t> const& embedded_surface_vertices() const;
    std::vector<meshless_sph_surface_vertex_t>& embedded_surface_vertices();

    meshless_sph_body_t* mechanical_model();
    meshless_sph_body_t const* mechanical_model() const;

  private:
    void prepare_vertices_for_surface_rendering();

    std::vector<vertex_type>
        render_vertices_; ///< The vertices_ member of tetrahedral_mesh_boundary_t will be
                          ///< considered as vertex positions in material space, while the
                          ///< render_vertices_ are the world space vertex positions
    std::vector<triangle_type> triangles_;

    std::vector<meshless_sph_surface_vertex_t> vertices_;
    meshless_sph_body_t* mechanical_model_;
};

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_MESHLESS_SURFACE_H