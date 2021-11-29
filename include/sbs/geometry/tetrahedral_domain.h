#ifndef SBS_GEOMETRY_TETRAHEDRAL_DOMAIN_H
#define SBS_GEOMETRY_TETRAHEDRAL_DOMAIN_H

#include "sbs/math/mapping.h"
#include "sbs/topology/tetrahedron_set.h"

#include <Discregrid/acceleration/bounding_sphere.hpp>
#include <Discregrid/acceleration/kd_tree.hpp>
#include <Eigen/Core>
#include <utility>
#include <vector>

namespace sbs {
namespace geometry {

scalar_type tetrahedron_volume(
    Eigen::Vector3d const& p1,
    Eigen::Vector3d const& p2,
    Eigen::Vector3d const& p3,
    Eigen::Vector3d const& p4);

class tetrahedral_domain_t
{
  public:
    class in_tetrahedron_query_t;

    // Constructors
    tetrahedral_domain_t() = default;
    tetrahedral_domain_t(
        std::vector<Eigen::Vector3d> const& points,
        std::vector<index_type> const& indices,
        scalar_type const query_error = 1e-2);

    tetrahedral_domain_t(tetrahedral_domain_t const& other);

    tetrahedral_domain_t& tetrahedral_domain_t::operator=(tetrahedral_domain_t const& other);

    // Accessors
    topology::tetrahedron_t const& tetrahedron(index_type ti) const;
    Eigen::Vector3d const& position(index_type i) const;
    topology::tetrahedron_set_t const& topology() const;
    index_type in_tetrahedron(Eigen::Vector3d const& X) const;
    bool contains(Eigen::Vector3d const& X) const;
    math::tetrahedron_barycentric_mapping_t const& barycentric_map(index_type ti) const;

  private:
    std::vector<Eigen::Vector3d> positions_; ///< Vertex positions
    topology::tetrahedron_set_t mesh_;       ///< Topology
    std::vector<math::tetrahedron_barycentric_mapping_t>
        tet_maps_; ///< Barycentric mapping per tetrahedron

  public:
    class in_tetrahedron_query_t : public Discregrid::KDTree<Discregrid::BoundingSphere>
    {
      public:
        using base_type = Discregrid::KDTree<Discregrid::BoundingSphere>;

        in_tetrahedron_query_t();
        in_tetrahedron_query_t(
            topology::tetrahedron_set_t const* topology,
            std::vector<Eigen::Vector3d> const* mesh_nodes,
            std::vector<math::tetrahedron_barycentric_mapping_t> const* tet_maps,
            scalar_type tolerance);

        in_tetrahedron_query_t(in_tetrahedron_query_t const& other) = default;
        in_tetrahedron_query_t(in_tetrahedron_query_t&& other)      = default;

        in_tetrahedron_query_t& operator=(in_tetrahedron_query_t const& other) = default;
        in_tetrahedron_query_t& operator=(in_tetrahedron_query_t&& other) noexcept = default;

        index_type in_tetrahedron(Eigen::Vector3d const& p) const;

        virtual Eigen::Vector3d entityPosition(unsigned int i) const override final;
        virtual void computeHull(unsigned int b, unsigned int n, Discregrid::BoundingSphere& hull)
            const override final;

        friend class tetrahedral_domain_t;

      private:
        topology::tetrahedron_set_t const* topology_;
        std::vector<Eigen::Vector3d> const* mesh_nodes_;
        std::vector<math::tetrahedron_barycentric_mapping_t> const* tet_maps_;
        std::vector<Eigen::Vector3d> tetrahedron_centers_;
        scalar_type tolerance_;
    };

  private:
    in_tetrahedron_query_t in_tetrahedron_query_; ///< Spatial searcher
};

std::pair<std::vector<Eigen::Vector3d>, std::vector<index_type>>
boundary_surface(tetrahedral_domain_t const& domain);

} // namespace geometry
} // namespace sbs

#endif // SBS_GEOMETRY_TETRAHEDRAL_DOMAIN_H