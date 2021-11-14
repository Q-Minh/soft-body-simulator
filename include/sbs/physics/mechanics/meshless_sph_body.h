#ifndef SBS_PHYSICS_MECHANICS_MESHLESS_BODY_H
#define SBS_PHYSICS_MECHANICS_MESHLESS_BODY_H

#include "sbs/physics/body/body.h"
#include "sbs/physics/collision/bvh_model.h"
#include "sbs/physics/mechanics/meshless_sph_surface.h"
#include "sbs/topology/tetrahedron_set.h"

#include <Discregrid/acceleration/bounding_sphere_hierarchy.hpp>
#include <Eigen/Core>
#include <sbs/aliases.h>
#include <vector>

namespace sbs {
namespace common {

struct geometry_t;

} // namespace common
namespace physics {
namespace mechanics {

class meshless_sph_node_t;

class meshless_sph_body_range_searcher_t : public Discregrid::KDTree<Discregrid::BoundingSphere>
{
  public:
    using base_type = Discregrid::KDTree<Discregrid::BoundingSphere>;

    meshless_sph_body_range_searcher_t();
    meshless_sph_body_range_searcher_t(std::vector<meshless_sph_node_t> const* nodes);

    meshless_sph_body_range_searcher_t(meshless_sph_body_range_searcher_t const& other) = default;
    meshless_sph_body_range_searcher_t(meshless_sph_body_range_searcher_t&& other)      = default;

    meshless_sph_body_range_searcher_t&
    operator=(meshless_sph_body_range_searcher_t const& other) = default;
    meshless_sph_body_range_searcher_t&
    operator=(meshless_sph_body_range_searcher_t&& other) noexcept = default;

    std::vector<index_type> neighbours_of(index_type const ni) const;
    std::vector<index_type> neighbours_of(Eigen::Vector3d const& p, scalar_type const h) const;
    std::vector<index_type> is_in_node_domains(Eigen::Vector3d const& p) const;

    virtual Eigen::Vector3d entityPosition(unsigned int i) const override final;
    virtual void computeHull(unsigned int b, unsigned int n, Discregrid::BoundingSphere& hull)
        const override final;

  private:
    std::vector<meshless_sph_node_t> const* nodes_;
};

class meshless_sph_body_t : public body::body_t
{
  public:
    /**
     * @brief
     * Constructs a meshless sph model by taking in an initial
     * tetrahedral geometry, extracting its boundary surface,
     * and then sampling meshless particles of the axis-aligned
     * grid englobing this boundary surface. The grid's resolution
     * is given by the resolution parameter. This determines the
     * sampling rate of the particles inside the boundary surface.
     *
     * @param simulation
     * @param id
     * @param geometry The initial tetrahedral geometry
     * @param h The support (or smoothing length) of meshless meshless_nodes as a multiplier of the
     * grid's cells' dimensions. For example, if the grid's resolution is 10x10x10 in a 10x10x10
     * domain. The grid's cells' dimensions will be 1x1x1. The smoothing length will be computed as
     * h*1. Increase h to include more neighbours.
     * @param resolution The resolution of the particle sampling grid
     */
    meshless_sph_body_t(
        xpbd::simulation_t& simulation,
        index_type id,
        common::geometry_t const& geometry,
        scalar_type const h,
        std::array<unsigned int, 3u> const& resolution);

    virtual visual_model_type const& visual_model() const override;
    virtual collision_model_type const& collision_model() const override;
    virtual visual_model_type& visual_model() override;
    virtual collision_model_type& collision_model() override;
    virtual void update_visual_model() override;
    virtual void update_collision_model() override;
    virtual void update_physical_model() override;
    void transform(Eigen::Affine3d const& affine);

    std::vector<meshless_sph_node_t> const& nodes() const;
    std::vector<meshless_sph_node_t>& nodes();
    meshless_sph_surface_t const& surface_mesh() const;
    meshless_sph_surface_t& surface_mesh();
    collision::point_bvh_model_t const& bvh() const;

    meshless_sph_body_range_searcher_t const& range_searcher() const;
    topology::tetrahedron_set_t const& topology() const;
    topology::tetrahedron_set_t& topology();
    scalar_type h() const;

    void initialize_physical_model();
    void initialize_visual_model();
    void initialize_collision_model();

  private:
    std::vector<meshless_sph_node_t> meshless_nodes_;
    meshless_sph_body_range_searcher_t
        material_space_range_query_; ///< Used for querying neighbours in material space
    topology::tetrahedron_set_t
        volumetric_topology_; ///< Only used for boundary surface extraction, the mechanical
                              ///< representation is still just a set of meshless meshless_nodes
    meshless_sph_surface_t visual_model_;
    collision::point_bvh_model_t collision_model_;
    scalar_type h_; ///< Support radius of meshless meshless_nodes
};

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_MESHLESS_BODY_H