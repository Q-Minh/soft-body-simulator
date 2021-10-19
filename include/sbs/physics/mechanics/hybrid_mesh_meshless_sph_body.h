#ifndef SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_SPH_BODY_H
#define SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_SPH_BODY_H

#include <Discregrid/acceleration/bounding_sphere.hpp>
#include <Discregrid/acceleration/kd_tree.hpp>
#include <Eigen/Core>
#include <sbs/physics/body.h>
#include <sbs/physics/collision/bvh_model.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_surface.h>
#include <sbs/physics/topology.h>

namespace sbs {
namespace common {
struct geometry_t;
} // namespace common

namespace physics {
namespace mechanics {

class hybrid_mesh_meshless_sph_node_t;

namespace detail {
namespace hybrid_mesh_meshless_sph {

class meshless_node_range_searcher_t : public Discregrid::KDTree<Discregrid::BoundingSphere>
{
  public:
    using base_type = Discregrid::KDTree<Discregrid::BoundingSphere>;

    meshless_node_range_searcher_t();
    meshless_node_range_searcher_t(std::vector<hybrid_mesh_meshless_sph_node_t> const* nodes);

    meshless_node_range_searcher_t(meshless_node_range_searcher_t const& other) = default;
    meshless_node_range_searcher_t(meshless_node_range_searcher_t&& other)      = default;

    meshless_node_range_searcher_t&
    operator=(meshless_node_range_searcher_t const& other) = default;
    meshless_node_range_searcher_t&
    operator=(meshless_node_range_searcher_t&& other) noexcept = default;

    std::vector<index_type> neighbours_of(index_type const ni) const;
    std::vector<index_type> neighbours_of(Eigen::Vector3d const& p, scalar_type const h) const;

    virtual Eigen::Vector3d entityPosition(unsigned int i) const override final;
    virtual void computeHull(unsigned int b, unsigned int n, Discregrid::BoundingSphere& hull)
        const override final;

  private:
    std::vector<hybrid_mesh_meshless_sph_node_t> const* nodes_;
};

class mesh_tetrahedron_range_searcher_t : public Discregrid::KDTree<Discregrid::BoundingSphere>
{
  public:
    using base_type = Discregrid::KDTree<Discregrid::BoundingSphere>;

    mesh_tetrahedron_range_searcher_t();
    mesh_tetrahedron_range_searcher_t(
        tetrahedron_set_t const* topology,
        std::vector<Eigen::Vector3d> const* mesh_nodes,
        std::vector<Eigen::Matrix4d> const* Ainv);

    mesh_tetrahedron_range_searcher_t(mesh_tetrahedron_range_searcher_t const& other) = default;
    mesh_tetrahedron_range_searcher_t(mesh_tetrahedron_range_searcher_t&& other)      = default;

    mesh_tetrahedron_range_searcher_t&
    operator=(mesh_tetrahedron_range_searcher_t const& other) = default;
    mesh_tetrahedron_range_searcher_t&
    operator=(mesh_tetrahedron_range_searcher_t&& other) noexcept = default;

    index_type in_tetrahedron(Eigen::Vector3d const& p) const;

    virtual Eigen::Vector3d entityPosition(unsigned int i) const override final;
    virtual void computeHull(unsigned int b, unsigned int n, Discregrid::BoundingSphere& hull)
        const override final;

  private:
    tetrahedron_set_t const* topology_;
    std::vector<Eigen::Vector3d> const* mesh_nodes_;
    std::vector<Eigen::Matrix4d> const* Ainv_;
    std::vector<Eigen::Vector3d> tetrahedron_centers_;
};

} // namespace hybrid_mesh_meshless_sph
} // namespace detail

class hybrid_mesh_meshless_sph_body_t : public physics::body_t
{
  public:
    /**
     * @brief
     * Constructs a hybrid mesh/meshless model. Takes all
     * boundary tets of the given tetrahedral mesh geometry,
     * and inserts particles there. The coupling is done by
     * having these boundary tets remove shape functions
     * of boundary mesh meshless_nodes. In addition, truncation of
     * the meshless meshless_nodes' domains is done to make sure
     * the meshless shape functions vanish at every interior
     * mesh node. This ensure continuity at the mixed boundary.
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
    hybrid_mesh_meshless_sph_body_t(
        simulation_t& simulation,
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
    virtual void transform(Eigen::Affine3d const& affine) override;

    std::vector<hybrid_mesh_meshless_sph_node_t> const& meshless_nodes() const;
    std::vector<hybrid_mesh_meshless_sph_node_t>& meshless_nodes();
    std::size_t mesh_node_count() const;
    std::size_t meshless_node_count() const;
    hybrid_mesh_meshless_sph_surface_t const& surface_mesh() const;
    hybrid_mesh_meshless_sph_surface_t& surface_mesh();
    collision::point_bvh_model_t const& bvh() const;

    std::size_t mixed_meshless_node_count() const;
    std::size_t interior_tetrahedron_count() const;
    std::size_t mesh_shape_function_count() const;

    tetrahedron_set_t const& topology() const;
    tetrahedron_set_t& topology();
    scalar_type h() const;

    /**
     *  I would like not to have to expose these range searchers that are in the
     * detail namespace, but for the moment, this will do.
     */

    detail::hybrid_mesh_meshless_sph::mesh_tetrahedron_range_searcher_t const&
    mesh_node_range_searcher() const;

    detail::hybrid_mesh_meshless_sph::meshless_node_range_searcher_t const&
    meshless_node_range_searcher() const;

    bool is_boundary_mesh_tetrahedron(index_type const ti) const;
    bool is_boundary_mesh_vertex(index_type const vi) const;
    Eigen::Matrix4d const& phi_i(index_type const ti) const;
    Eigen::Matrix<scalar_type, 4, 3> grad_phi_i(index_type const ti) const;
    Eigen::Vector3d grad_phi_i(index_type const ti, std::uint8_t v) const;

    void initialize_physical_model();
    void initialize_visual_model();
    void initialize_collision_model();

  private:
    std::vector<Eigen::Vector3d>
        mesh_nodes_; ///< Material space positions of the initial tetrahedral geometry
    std::vector<Eigen::Vector3d>
        world_space_positions_; ///< Positions of each mesh node in world space
    std::vector<hybrid_mesh_meshless_sph_node_t> meshless_nodes_;
    detail::hybrid_mesh_meshless_sph::meshless_node_range_searcher_t
        material_space_meshless_node_searcher_; ///< Used for querying neighbours in material space
    detail::hybrid_mesh_meshless_sph::mesh_tetrahedron_range_searcher_t
        material_space_mesh_node_searcher_;
    tetrahedron_set_t volumetric_topology_; ///< Is not the topology of the mechanical model.
    std::vector<bool>
        is_boundary_vertex_; ///< vi is a boundary vertex if is_boundary_vertex_ == true
    std::vector<bool>
        is_boundary_tetrahedron_; ///< ti is a boundary tet if is_boundary_tetrahedron_[ti] == true
    std::vector<Eigen::Matrix4d>
        Ainv_; ///< the 4 barycentric coordinates over a given tetrahedron is given by Ainv * X.
               ///< Each row of the resulting vector of barycentric coordinates corresponds to the
               ///< shape function associated with vertex vi \in 0,1,2,3 in ti.

    hybrid_mesh_meshless_sph_surface_t visual_model_;
    collision::point_bvh_model_t collision_model_;
    scalar_type h_; ///< Support radius of meshless meshless_nodes
};

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_SPH_BODY_H