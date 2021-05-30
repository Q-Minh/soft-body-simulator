#ifndef SBS_PHYSICS_CUTTING_CUT_TETRAHEDRON_H
#define SBS_PHYSICS_CUTTING_CUT_TETRAHEDRON_H

#include "physics/collision/intersections.h"

#include <Eigen/Core>
#include <array>
#include <utility>

namespace sbs {
namespace physics {
namespace cutting {

class tetrahedron_mesh_cutter_t
{
  public:
    using tetrahedron_type      = Eigen::Matrix<std::uint32_t, 4, 1>;
    using tetrahedra_type       = Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic>;
    using positions_type        = Eigen::Matrix3Xd;
    using tetrahedral_mesh_type = std::pair<positions_type, tetrahedra_type>;

    std::vector<tetrahedral_mesh_type> subdivide_mesh(
        std::byte const& edge_intersection_mask,
        positions_type const& TV,
        tetrahedra_type const& TT,
        std::uint32_t tetrahedron,
        std::array<Eigen::Vector3d, 6u> const& edge_intersection_points,
        std::array<Eigen::Vector3d, 4u> const& face_intersection_points);

  private:
    std::vector<tetrahedral_mesh_type> subdivide_mesh_for_common_case_1(
        positions_type const& TV,
        tetrahedra_type const& TT,
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<Eigen::Vector3d, 3u> const& edge_intersection_points);

    std::vector<tetrahedral_mesh_type> subdivide_mesh_for_common_case_2(
        positions_type const& TV,
        tetrahedra_type const& TT,
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<Eigen::Vector3d, 4u> const& edge_intersection_points);

    std::vector<tetrahedral_mesh_type> subdivide_mesh_for_common_case_3(
        positions_type const& TV,
        tetrahedra_type const& TT,
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<Eigen::Vector3d, 1u> const& edge_intersection_points,
        std::array<Eigen::Vector3d, 2u> const& face_intersection_points);

    std::vector<tetrahedral_mesh_type> subdivide_mesh_for_common_case_4(
        positions_type const& TV,
        tetrahedra_type const& TT,
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<Eigen::Vector3d, 2u> const& edge_intersection_points,
        std::array<Eigen::Vector3d, 2u> const& face_intersection_points);

    std::vector<tetrahedral_mesh_type> subdivide_mesh_for_common_case_5(
        positions_type const& TV,
        tetrahedra_type const& TT,
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<Eigen::Vector3d, 3u> const& edge_intersection_points,
        std::array<Eigen::Vector3d, 2u> const& face_intersection_points,
        bool symmetry = false);
};

std::vector<std::pair<Eigen::Matrix3Xd, Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic>>>
cut_tetrahedron(
    Eigen::Matrix3Xd const& V,
    Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic> const& T,
    std::uint32_t tetrahedron,
    std::byte const edge_intersection_mask,
    std::array<Eigen::Vector3d, 6u> const& edge_intersections,
    std::array<Eigen::Vector3d, 4u> const& face_intersections);

} // namespace cutting
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_CUTTING_CUT_TETRAHEDRON_H