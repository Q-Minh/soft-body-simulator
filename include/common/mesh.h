#ifndef SBS_COMMON_MESH_H
#define SBS_COMMON_MESH_H

#include "geometry.h"

#include <Eigen/Core>
#include <vector>

namespace sbs {
namespace common {

/**
 * @brief Tetrahedral mesh type
 */
class shared_vertex_mesh_t
{
  public:
    using index_type       = std::uint32_t;
    using positions_type   = Eigen::Matrix3Xd;
    using elements_type    = Eigen::Matrix<index_type, Eigen::Dynamic, Eigen::Dynamic>;
    using triangles_type   = Eigen::Matrix<index_type, 3, Eigen::Dynamic>;
    using masses_type      = Eigen::VectorXd;
    using velocities_type  = Eigen::Matrix3Xd;
    using forces_type      = Eigen::Matrix3Xd;
    using tetrahedron_type = Eigen::Matrix<index_type, 4, 1>;
    using tetrahedra_type  = Eigen::Matrix<index_type, 4, Eigen::Dynamic>;
    using triangle_type    = Eigen::Matrix<index_type, 3, 1>;

    using uv_coordinates_type = Eigen::Matrix2Xf;
    using normals_type        = Eigen::Matrix3Xd;
    using colors_type         = Eigen::Matrix3Xf;

    shared_vertex_mesh_t() = default;
    shared_vertex_mesh_t(common::geometry_t const& geometry);
    shared_vertex_mesh_t(positions_type const& P, tetrahedra_type const& T);

    void rescale(
        Eigen::Vector3d const& boxmin = Eigen::Vector3d{-1., -1., -1.},
        Eigen::Vector3d const& boxmax = Eigen::Vector3d{+1., +1., +1.});

    positions_type const& positions() const;
    positions_type& positions();

    elements_type const& elements() const;
    elements_type& elements();

    triangles_type const& faces() const;
    triangles_type& faces();

    masses_type const& masses() const;
    masses_type& masses();

    velocities_type const& velocities() const;
    velocities_type& velocities();

    forces_type const& forces() const;
    forces_type& forces();

    /**
     * Rendering
     */
    positions_type const& vertices() const;
    positions_type& vertices();

    uv_coordinates_type const& uvs() const;
    uv_coordinates_type& uvs();

    normals_type const& normals() const;
    normals_type& normals();

    colors_type const& colors() const;
    colors_type& colors();

    void set_color(Eigen::Vector3f const rgb);

    void extract_boundary_surface_mesh();
    void extract_boundary_normals();

  private:
    positions_type positions_;
    elements_type elements_;
    masses_type masses_;
    velocities_type velocities_;
    forces_type forces_;
    uv_coordinates_type uvs_;
    colors_type colors_;

    triangles_type boundary_faces_;
    positions_type boundary_vertices_;
    uv_coordinates_type boundary_uvs_;
    colors_type boundary_colors_;
    normals_type normals_;
};

std::vector<std::pair<std::uint32_t, std::uint32_t>> edges(shared_vertex_mesh_t const& mesh);

} // namespace common
} // namespace sbs

#endif // SBS_COMMON_MESH_H