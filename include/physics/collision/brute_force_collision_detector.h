#ifndef SBS_PHYSICS_COLLISION_BRUTE_FORCE_COLLISION_DETECTOR_H
#define SBS_PHYSICS_COLLISION_BRUTE_FORCE_COLLISION_DETECTOR_H

#include "physics/collision/collision_detector.h"
#include "physics/mesh.h"

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
    : public collision_detector_i<std::pair<SimulatedMeshType*, std::uint32_t>>
{
  public:
    using mesh_type      = SimulatedMeshType;
    using base_mesh_type = physics::simulated_mesh_i<mesh_type>;
    using element_type   = std::pair<mesh_type* /*body*/, std::uint32_t /*ti*/>;

    brute_force_collision_detector_t() = default;

    template <class BaseMeshTypePointerIt>
    explicit brute_force_collision_detector_t(
        BaseMeshTypePointerIt first,
        BaseMeshTypePointerIt last);

    void add_mesh(base_mesh_type* mesh);

    virtual std::list<element_type> intersect(common::tetrahedron_t const& t) const override;
    virtual std::list<element_type> intersect(common::triangle_t const& t) const override;
    virtual element_type intersect(common::point_t const& p) const override;

  private:
    template <class GeometricPrimitiveType>
    std::list<element_type> intersect_base(GeometricPrimitiveType const& geometric_primitive) const;

  private:
    std::vector<base_mesh_type*> bodies_;
};

template <class SimulatedMeshType>
template <class BaseMeshTypePointerIt>
inline brute_force_collision_detector_t<SimulatedMeshType>::brute_force_collision_detector_t(
    BaseMeshTypePointerIt first,
    BaseMeshTypePointerIt last)
    : bodies_{first, last}
{
}

template <class SimulatedMeshType>
template <class GeometricPrimitiveType>
inline std::list<brute_force_collision_detector_t<SimulatedMeshType>::element_type>
brute_force_collision_detector_t<SimulatedMeshType>::intersect_base(
    GeometricPrimitiveType const& geometric_primitive) const
{
    std::vector<std::list<element_type>> intersections(bodies_.size(), std::list<element_type>{});

    auto const compute_intersections_for_body = [geometric_primitive](base_mesh_type* body) {
        std::list<element_type> intersected_elements{};

        std::vector<tetrahedron_t> const& tetrahedra   = body->tetrahedra();
        std::vector<physics::vertex_t> const& vertices = body->vertices();

        for (std::size_t ti = 0u; i < tetrahedra.size(); ++ti)
        {
            tetrahedron_t const& tetrahedron = tetrahedra[ti];

            index_type const v1 = tetrahedron.v1();
            index_type const v2 = tetrahedron.v2();
            index_type const v3 = tetrahedron.v3();
            index_type const v4 = tetrahedron.v4();

            auto const& p1 = vertices[v1].position();
            auto const& p2 = vertices[v2].position();
            auto const& p3 = vertices[v3].position();
            auto const& p4 = vertices[v4].position();

            common::tetrahedron_t const tetrahedron_primitive{p1, p2, p3, p4};
            if (common::intersects(geometric_primitive, tetrahedron_primitive))
            {
                intersected_elements.push_back({static_cast<mesh_type*>(body), ti});
            }
        }
        return intersected_elements;
    };

    std::transform(
        bodies_.begin(),
        bodies_.end(),
        intersections.begin(),
        compute_intersections_for_body);

    auto const reduce_op = [](std::list<element_type>&& accum, std::list<element_type>& next) {
        accum.splice(accum.begin(), next);
        return std::move(accum);
    };
    auto const all_intersection_pairs = std::reduce(
        intersections.begin(),
        intersections.end(),
        std::list<element_type>{},
        reduce_op);
    return all_intersection_pairs;
}

template <class SimulatedMeshType>
inline void brute_force_collision_detector_t<SimulatedMeshType>::add_mesh(base_mesh_type* mesh)
{
    bodies_.push_back(mesh);
}

template <class SimulatedMeshType>
inline std::list<brute_force_collision_detector_t<SimulatedMeshType>::element_type>
brute_force_collision_detector_t<SimulatedMeshType>::intersect(common::tetrahedron_t const& t) const
{
    return intersect_base<common::tetrahedron_t>(t);
}

template <class SimulatedMeshType>
inline std::list<brute_force_collision_detector_t<SimulatedMeshType>::element_type>
brute_force_collision_detector_t<SimulatedMeshType>::intersect(common::triangle_t const& t) const
{
    return intersect_base<common::triangle_t>(t);
}

template <class SimulatedMeshType>
inline brute_force_collision_detector_t<SimulatedMeshType>::element_type
brute_force_collision_detector_t<SimulatedMeshType>::intersect(common::point_t const& p) const
{
    return intersect_base<common::point_t>(p).front();
}

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_BRUTE_FORCE_COLLISION_DETECTOR_H