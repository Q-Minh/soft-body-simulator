#ifndef SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_SPH_SURFACE_H
#define SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_SPH_SURFACE_H

#include <Eigen/Core>
#include <sbs/physics/tetrahedral_mesh_boundary.h>

namespace sbs {
namespace physics {
namespace mechanics {

class hybrid_mesh_meshless_sph_body_t;

class hybrid_mesh_meshless_sph_surface_t : public tetrahedral_mesh_boundary_t
{
  public:
    using vertex_type   = tetrahedral_mesh_boundary_t::vertex_type;
    using triangle_type = tetrahedral_mesh_boundary_t::triangle_type;

    hybrid_mesh_meshless_sph_surface_t() = default;
    hybrid_mesh_meshless_sph_surface_t(hybrid_mesh_meshless_sph_body_t* mechanical_model);

    hybrid_mesh_meshless_sph_surface_t(hybrid_mesh_meshless_sph_surface_t const& other) = default;
    hybrid_mesh_meshless_sph_surface_t(hybrid_mesh_meshless_sph_surface_t&& other)      = default;
    hybrid_mesh_meshless_sph_surface_t&
    operator=(hybrid_mesh_meshless_sph_surface_t const& other) = default;
    hybrid_mesh_meshless_sph_surface_t&
    operator=(hybrid_mesh_meshless_sph_surface_t&& other) = default;

    void initialize_interpolation_scheme();

    virtual vertex_type vertex(std::size_t vi) const override;
    virtual void prepare_vertices_for_rendering() override;

    vertex_type& world_space_vertex(std::size_t vi);
    vertex_type& material_space_vertex(std::size_t vi);
    void compute_positions();
    void compute_normals();

  private:
    void prepare_vertices_for_surface_rendering();

    std::vector<vertex_type>
        world_space_vertices_; ///< The vertices_ member of tetrahedral_mesh_boundary_t will be
                               ///< considered as vertex positions in material space, while the
                               ///< world_space_vertices_ are the world space vertex positions

    std::vector<std::vector<Eigen::Vector3d>> Xkjs_;
    std::vector<std::vector<scalar_type>> Wkjs_;
    std::vector<scalar_type> sks_;
    std::vector<std::vector<index_type>> neighbours_;
    hybrid_mesh_meshless_sph_body_t* mechanical_model_;
};

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_SPH_SURFACE_H