#ifndef SBS_PHYSICS_CUTTING_CUT_TETRAHEDRON_H
#define SBS_PHYSICS_CUTTING_CUT_TETRAHEDRON_H

#include "common/mesh.h"
#include "physics/collision/intersections.h"

#include <Eigen/Core>
#include <array>
#include <list>
#include <map>
#include <optional>
#include <utility>

namespace sbs {
namespace physics {
namespace cutting {

class tetrahedral_mesh_cutter_t
{
  public:
    using tetrahedron_type        = Eigen::Matrix<std::uint32_t, 4, 1>;
    using tetrahedra_type         = Eigen::Matrix<std::uint32_t, 4, Eigen::Dynamic>;
    using positions_type          = Eigen::Matrix3Xd;
    using masses_type             = Eigen::VectorXd;
    using velocities_type         = Eigen::Matrix3Xd;
    using forces_type             = Eigen::Matrix3Xd;
    using subdivided_element_type = tetrahedra_type;

    /**
     * @brief
     * @param mesh
     */
    tetrahedral_mesh_cutter_t(common::shared_vertex_mesh_t& mesh);
    tetrahedral_mesh_cutter_t(tetrahedral_mesh_cutter_t const& other)     = default;
    tetrahedral_mesh_cutter_t(tetrahedral_mesh_cutter_t&& other) noexcept = default;

    /**
     * @brief
     * @param cutting_surface
     * @return
     */
    bool cut_tetrahedral_mesh(common::shared_vertex_surface_mesh_t const& cutting_surface);

    /**
     * @brief
     */
    void finalize_cut();

  protected:
    using triangle_facet_key_type = std::array<std::uint32_t, 3u>;
    using edge_facet_key_type     = std::array<std::uint32_t, 2u>;
    using tetrahedron_key_type    = std::uint32_t;
    using cutting_mask_type       = std::byte;

    struct cut_facet_t
    {
        std::uint32_t vi{0u}; ///< vi is the index of the newly created vertex that will be used for
                              ///< subdivided elements under the cutting surface, while vi + 1u is
                              ///< the index of the newly created vertex that will be used
                              ///< for the subdivided elements over the surface
        bool has_created_new_vertex{false};
        std::uint32_t v1, v2; ///< represents the cut edge (v1,v2) where v1 is under the cutting
                              ///< surface, v2 is over the cutting surface
        Eigen::Vector3d intersection{};
    };

    bool intersects(
        std::uint32_t tetrahedron,
        common::shared_vertex_surface_mesh_t const& cutting_surface) const;

  private:
    struct edge_facet_key_less
    {
        bool operator()(edge_facet_key_type const& e1, edge_facet_key_type const& e2) const
        {
            edge_facet_key_type e1p{e1};
            edge_facet_key_type e2p{e2};
            std::sort(e1p.begin(), e1p.end());
            std::sort(e2p.begin(), e2p.end());
            return e1p < e2p;
        }
    };

    struct triangle_facet_key_less
    {
        bool operator()(triangle_facet_key_type const& t1, triangle_facet_key_type const& t2) const
        {
            triangle_facet_key_type t1p{t1};
            triangle_facet_key_type t2p{t2};
            std::sort(t1p.begin(), t1p.end());
            std::sort(t2p.begin(), t2p.end());
            return t1p < t2p;
        }
    };

    void detect_cuts(common::shared_vertex_surface_mesh_t const& cutting_surface);
    void detect_cuts(
        std::uint32_t tetrahedron,
        common::shared_vertex_surface_mesh_t const& cutting_surface);

    /**
     * @brief
     * Create new vertices along with their physical quantities (mass, force, velocity)
     */
    void create_new_vertices();

    /**
     * @brief
     * Perform mesh subdivision based on current and past intersections
     * @param cutting_surface
     * @return
     */
    bool update_topology(common::shared_vertex_surface_mesh_t const& cutting_surface);

    std::optional<subdivided_element_type>
    subdivide_mesh(std::byte const& edge_intersection_mask, std::uint32_t tetrahedron);

    subdivided_element_type subdivide_mesh_for_common_case_1(
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<edge_facet_key_type, 3u> const& edge_intersections);

    subdivided_element_type subdivide_mesh_for_common_case_2(
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<edge_facet_key_type, 4u> const& edge_intersections);

    subdivided_element_type subdivide_mesh_for_common_case_3(
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<edge_facet_key_type, 1u> const& edge_intersections,
        std::array<triangle_facet_key_type, 2u> const& face_intersections);

    subdivided_element_type subdivide_mesh_for_common_case_4(
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<edge_facet_key_type, 2u> const& edge_intersections,
        std::array<triangle_facet_key_type, 2u> const& face_intersections);

    subdivided_element_type subdivide_mesh_for_common_case_5(
        std::uint32_t tetrahedron,
        std::array<std::uint32_t, 4u> const& vertex_ordering,
        std::array<edge_facet_key_type, 3u> const& edge_intersections,
        std::array<triangle_facet_key_type, 2u> const& face_intersections,
        bool symmetry = false);

    common::shared_vertex_mesh_t& mesh_;
    std::map<edge_facet_key_type, cut_facet_t, edge_facet_key_less> cut_edges_;
    std::map<triangle_facet_key_type, cut_facet_t, triangle_facet_key_less> cut_faces_;
    std::map<tetrahedron_key_type, cutting_mask_type> cut_tetrahedra_;
    std::uint32_t previous_vertex_count_;
};

} // namespace cutting
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_CUTTING_CUT_TETRAHEDRON_H