#ifndef SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H

#include "constraint.h"

#include <Eigen/Core>

namespace sbs {
namespace physics {
namespace xpbd {

class green_constraint_t : public constraint_t
{
  public:
    using scalar_type   = typename constraint_t::scalar_type;
    using index_type    = std::uint32_t;
    using body_ptr_type = tetrahedral_mesh_t*;

    green_constraint_t(
        scalar_type const alpha,
        body_ptr_type b,
        index_type ti,
        scalar_type young_modulus,
        scalar_type poisson_ratio);

    virtual void project(
        std::vector<std::shared_ptr<tetrahedral_mesh_t>> const& bodies,
        scalar_type& lagrange_multiplier,
        scalar_type const dt) const override;

  protected:
    scalar_type signed_volume(
        Eigen::Vector3d const& p1,
        Eigen::Vector3d const& p2,
        Eigen::Vector3d const& p3,
        Eigen::Vector3d const& p4) const;

  private:
    body_ptr_type b_;
    index_type ti_;

    Eigen::Matrix3d DmInv_;
    scalar_type V0_;
    scalar_type mu_;
    scalar_type lambda_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H