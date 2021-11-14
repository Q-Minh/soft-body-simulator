#ifndef SBS_PHYSICS_XPBD_HYBRID_MESH_MESHLESS_MLS_STVK_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_HYBRID_MESH_MESHLESS_MLS_STVK_CONSTRAINT_H

#include "sbs/physics/xpbd/constraint.h"

#include <Eigen/Core>
#include <optional>

namespace sbs {
namespace physics {
namespace mechanics {

class hybrid_mesh_meshless_mls_node_t;

} // namespace mechanics

namespace xpbd {

class hybrid_mesh_meshless_mls_constraint_t : public constraint_t
{
  public:
    hybrid_mesh_meshless_mls_constraint_t(
        scalar_type const alpha,
        scalar_type const beta,
        simulation_t const& simulation,
        index_type bi,
        index_type ni,
        scalar_type young_modulus,
        scalar_type poisson_ratio,
        mechanics::hybrid_mesh_meshless_mls_node_t& node);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

    Eigen::Matrix3d deformation_gradient(simulation_t& simulation) const;
    Eigen::Matrix3d green_strain(Eigen::Matrix3d const& Fi) const;
    std::pair<scalar_type, Eigen::Matrix3d>
    strain_energy_and_stress(Eigen::Matrix3d const& Fi, Eigen::Matrix3d const& Ei) const;
    scalar_type const C(scalar_type const Psi) const;
    std::vector<Eigen::Vector3d> dCdxk(Eigen::Matrix3d const& dPsidFi) const;

    std::array<std::optional<Eigen::Vector3d /*gradient w.r.t. mesh node i*/>, 4u>
    dCdvi(Eigen::Matrix3d const& dPsidFi) const;

  private:
    mechanics::hybrid_mesh_meshless_mls_node_t& node_;
    index_type bi_;
    index_type ni_;
    scalar_type mu_;
    scalar_type lambda_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_HYBRID_MESH_MESHLESS_MLS_STVK_CONSTRAINT_H