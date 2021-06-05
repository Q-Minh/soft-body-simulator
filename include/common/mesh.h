#ifndef SBS_COMMON_MESH_H
#define SBS_COMMON_MESH_H

#include "geometry.h"

#include <Eigen/Core>
#include <vector>

namespace sbs {
namespace common {

/**
 * @brief Triangle mesh type
 */
class shared_vertex_surface_mesh_t
{
  public:
    using index_type     = std::uint32_t;
    using vertices_type  = Eigen::Matrix3Xd;
    using triangles_type = Eigen::Matrix<index_type, 3, Eigen::Dynamic>;
    using triangle_type  = Eigen::Matrix<index_type, 3, 1>;
    using normals_type   = Eigen::Matrix3Xd;
    using colors_type    = Eigen::Matrix3Xf;
    using index_map_type = std::vector<std::uint32_t>;
    using edge_type      = std::pair<index_type, index_type>;

    shared_vertex_surface_mesh_t() = default;
    shared_vertex_surface_mesh_t(common::geometry_t const& geometry);

    void set_color(Eigen::Vector3f const rgb = Eigen::Vector3f{1.f, 1.f, 0.f});
    void extract_normals();

    std::vector<edge_type> boundary_edges() const;

    vertices_type const& vertices() const;
    vertices_type& vertices();

    normals_type const& normals() const;
    normals_type& normals();

    colors_type const& colors() const;
    colors_type& colors();

    triangles_type const& triangles() const;
    triangles_type& triangles();

    index_map_type const& index_map() const;
    index_map_type& index_map();

  private:
    vertices_type vertices_;
    normals_type normals_;
    colors_type colors_;
    triangles_type triangles_;
    index_map_type
        index_map_; ///< Maps this surface mesh's vertex indices to their original
                    ///< tetrahedral mesh's vertex indices. This is useful to relate
                    ///< this surface mesh as a boundary to an original tetrahedral mesh,
                    ///< while not having to send all of the original tetrahedral mesh's vertices
                    ///< to the GPU for rendering. By having the index map, this surface mesh
                    ///< does not need to store all of the original tetrahedral mesh's vertices.
};

/**
 * @brief Tetrahedral mesh type for soft body physics (includes forces, velocities and masses)
 */
class shared_vertex_mesh_t
{
  public:
    using index_type       = std::uint32_t;
    using positions_type   = Eigen::Matrix3Xd;
    using tetrahedra_type  = Eigen::Matrix<index_type, 4, Eigen::Dynamic>;
    using tetrahedron_type = Eigen::Matrix<index_type, 4, 1>;
    using elements_type    = tetrahedra_type;
    using masses_type      = Eigen::VectorXd;
    using velocities_type  = Eigen::Matrix3Xd;
    using forces_type      = Eigen::Matrix3Xd;

    shared_vertex_mesh_t() = default;
    shared_vertex_mesh_t(common::geometry_t const& geometry);
    shared_vertex_mesh_t(positions_type const& P, tetrahedra_type const& T);
    shared_vertex_mesh_t(
        positions_type const& P,
        tetrahedra_type const& T,
        masses_type const& M,
        velocities_type const& V,
        forces_type const& F);

    positions_type const& positions() const;
    positions_type& positions();

    elements_type const& elements() const;
    elements_type& elements();

    masses_type const& masses() const;
    masses_type& masses();

    velocities_type const& velocities() const;
    velocities_type& velocities();

    forces_type const& forces() const;
    forces_type& forces();

    shared_vertex_surface_mesh_t
    boundary_surface_mesh(Eigen::Vector3f const& color = Eigen::Vector3f{1.f, 1.f, 0.f}) const;

    shared_vertex_surface_mesh_t
    facets(Eigen::Vector3f const& color = Eigen::Vector3f{1.f, 1.f, 0.f}) const;

  private:
    positions_type positions_;
    elements_type elements_;
    masses_type masses_;
    velocities_type velocities_;
    forces_type forces_;
};

std::vector<std::pair<std::uint32_t, std::uint32_t>> edges(shared_vertex_mesh_t const& mesh);

Eigen::Matrix3Xd rescale(
    Eigen::Matrix3Xd const& positions,
    Eigen::Vector3d const& boxmin = Eigen::Vector3d{-1., -1., -1.},
    Eigen::Vector3d const& boxmax = Eigen::Vector3d{+1., +1., +1.});

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_MESH_H