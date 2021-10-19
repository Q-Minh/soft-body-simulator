#ifndef SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_SPH_NODE_H
#define SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_SPH_NODE_H

#include <sbs/physics/mechanics/meshless_sph_node.h>

namespace sbs {
namespace physics {
namespace mechanics {

class hybrid_mesh_meshless_sph_body_t;

class hybrid_mesh_meshless_sph_node_t : public meshless_sph_node_t
{
  public:
    hybrid_mesh_meshless_sph_node_t(
        index_type const i,
        functions::poly6_kernel_t const& kernel,
        hybrid_mesh_meshless_sph_body_t& body);

    void initialize(
        Eigen::Vector3d const& xi,
        std::vector<Eigen::Vector3d const*> const& pj,
        std::vector<index_type> const& neighbours,
        index_type const ti);

    bool is_mixed_particle() const;

    hybrid_mesh_meshless_sph_body_t const& body() const;
    index_type ti() const;

  private:
    hybrid_mesh_meshless_sph_body_t& body_;
    index_type ti_; ///< Boundary tetrahedron in which this meshless node lives
};

} // namespace mechanics
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_MECHANICS_HYBRID_MESH_MESHLESS_SPH_NODE_H