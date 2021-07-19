#ifndef SBS_PHYSICS_COLLISION_BRUTE_FORCE_COLLISION_DETECTOR_H
#define SBS_PHYSICS_COLLISION_BRUTE_FORCE_COLLISION_DETECTOR_H

//#include "sbs/physics/collision/collision_detector.h"

#include "sbs/common/mesh.h"
#include "sbs/common/primitive.h"
#include "sbs/physics/mesh.h"

#include <list>
#include <numeric>
#include <vector>

namespace sbs {
namespace physics {

/**
 * @brief
 * @tparam ElementType Type of the elements we want to store.
 * @tparam SimulatedMeshType Type of the bodies to include in the collision detection system
 */
template <class SimulatedMeshType>
class brute_force_collision_detector_t
{
  public:
    using simulated_mesh_type    = SimulatedMeshType;
    using surface_mesh_type      = common::shared_vertex_surface_mesh_i;
    using intersection_pair_type = std::pair<std::uint32_t /*ti*/, common::triangle_t>;

    void add_environment_body(surface_mesh_type const* mesh);
    std::list<intersection_pair_type> intersect(simulated_mesh_type const* mesh) const;

  private:
    std::list<intersection_pair_type>
    intersect_single_triangle(simulated_mesh_type const* mesh, common::triangle_t const& t) const;

  private:
    std::vector<surface_mesh_type const*> environment_bodies_;
};

template <class SimulatedMeshType>
inline void brute_force_collision_detector_t<SimulatedMeshType>::add_environment_body(
    surface_mesh_type const* mesh)
{
    environment_bodies_.push_back(mesh);
}

template <class SimulatedMeshType>
inline std::list<
    typename brute_force_collision_detector_t<SimulatedMeshType>::intersection_pair_type>
brute_force_collision_detector_t<SimulatedMeshType>::intersect(
    simulated_mesh_type const* mesh) const
{
    std::list<intersection_pair_type> intersection_pairs{};
    for (auto const& env_body : environment_bodies_)
    {
        for (std::size_t fi = 0u; fi < env_body->triangle_count(); ++fi)
        {
            auto const f  = env_body->triangle(fi);
            auto const v1 = env_body->vertex(f.v1);
            auto const v2 = env_body->vertex(f.v2);
            auto const v3 = env_body->vertex(f.v3);

            common::point_t const p1{v1.x, v1.y, v1.z};
            common::point_t const p2{v2.x, v2.y, v2.z};
            common::point_t const p3{v3.x, v3.y, v3.z};

            common::triangle_t const triangle{p1, p2, p3};
            std::list<intersection_pair_type> triangle_tet_intersections =
                intersect_single_triangle(mesh, triangle);
            intersection_pairs.splice(intersection_pairs.begin(), triangle_tet_intersections);
        }
    }
    return intersection_pairs;
}

template <class SimulatedMeshType>
inline std::list<
    typename brute_force_collision_detector_t<SimulatedMeshType>::intersection_pair_type>
brute_force_collision_detector_t<SimulatedMeshType>::intersect_single_triangle(
    simulated_mesh_type const* mesh,
    common::triangle_t const& triangle) const
{
    std::list<intersection_pair_type> intersection_pairs{};
    std::vector<physics::tetrahedron_t> const& tetrahedra = mesh->tetrahedra();
    std::vector<physics::vertex_t> const& vertices        = mesh->vertices();
    for (std::size_t ti = 0u; ti < tetrahedra.size(); ++ti)
    {
        physics::tetrahedron_t const& t = tetrahedra[ti];

        auto const p1 = vertices[t.v1()].position();
        auto const p2 = vertices[t.v2()].position();
        auto const p3 = vertices[t.v3()].position();
        auto const p4 = vertices[t.v4()].position();

        common::tetrahedron_t const tetrahedron{p1, p2, p3, p4};
        if (common::intersects(triangle, tetrahedron))
        {
            intersection_pairs.push_back({static_cast<index_type>(ti), triangle});
        }
    }
    return intersection_pairs;
}

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_BRUTE_FORCE_COLLISION_DETECTOR_H