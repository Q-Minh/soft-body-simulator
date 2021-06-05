#ifndef SBS_PHYSICS_CUTTING_CUT_TETRAHEDRON_H
#define SBS_PHYSICS_CUTTING_CUT_TETRAHEDRON_H

#include "common/mesh.h"
#include "physics/collision/intersections.h"

#include <Eigen/Core>
#include <array>
#include <list>
#include <optional>
#include <utility>

namespace sbs {
namespace physics {
namespace cutting {

class tetrahedron_mesh_cutter_t
{
  public:
    using tetrahedron_type = Eigen::Matrix<std::uint32_t, 4, 1>;
    using tetrahedra_type  = Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic>;
    using positions_type   = Eigen::Matrix3Xd;
    using masses_type      = Eigen::VectorXd;
    using velocities_type  = Eigen::Matrix3Xd;
    using forces_type      = Eigen::Matrix3Xd;
    using subdivided_element_type =
        std::tuple<positions_type, tetrahedra_type, masses_type, velocities_type, forces_type>;

    std::optional<subdivided_element_type> subdivide_mesh(
        std::byte const& edge_intersection_mask,
        common::shared_vertex_mesh_t const& mesh,
        std::uint32_t tetrahedron,
        std::array<Eigen::Vector3d, 6u> const& edge_intersection_points,
        std::array<Eigen::Vector3d, 4u> const& face_intersection_points);

  private:
    subdivided_element_type subdivide_mesh_for_common_case_1(
        common::shared_vertex_mesh_t const& mesh,
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<Eigen::Vector3d, 3u> const& edge_intersection_points);

    subdivided_element_type subdivide_mesh_for_common_case_2(
        common::shared_vertex_mesh_t const& mesh,
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<Eigen::Vector3d, 4u> const& edge_intersection_points);

    subdivided_element_type subdivide_mesh_for_common_case_3(
        common::shared_vertex_mesh_t const& mesh,
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<Eigen::Vector3d, 1u> const& edge_intersection_points,
        std::array<Eigen::Vector3d, 2u> const& face_intersection_points);

    subdivided_element_type subdivide_mesh_for_common_case_4(
        common::shared_vertex_mesh_t const& mesh,
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<Eigen::Vector3d, 2u> const& edge_intersection_points,
        std::array<Eigen::Vector3d, 2u> const& face_intersection_points);

    subdivided_element_type subdivide_mesh_for_common_case_5(
        common::shared_vertex_mesh_t const& mesh,
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<Eigen::Vector3d, 3u> const& edge_intersection_points,
        std::array<Eigen::Vector3d, 2u> const& face_intersection_points,
        bool symmetry = false);
};

std::optional<tetrahedron_mesh_cutter_t::subdivided_element_type> cut_tetrahedron(
    common::shared_vertex_mesh_t const& mesh,
    std::uint32_t tetrahedron,
    std::byte const edge_intersection_mask,
    std::array<Eigen::Vector3d, 6u> const& edge_intersections,
    std::array<Eigen::Vector3d, 4u> const& face_intersections);

std::optional<tetrahedron_mesh_cutter_t::subdivided_element_type> cut_tetrahedron(
    common::shared_vertex_mesh_t const& mesh,
    std::uint32_t tetrahedron,
    common::shared_vertex_surface_mesh_t const& cutting_surface);

bool cut_tetrahedral_mesh(
    common::shared_vertex_mesh_t& mesh,
    common::shared_vertex_surface_mesh_t const& cutting_surface);

} // namespace cutting
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_CUTTING_CUT_TETRAHEDRON_H