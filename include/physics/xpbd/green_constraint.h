#ifndef SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H

#include "constraint.h"

#include <Eigen/Core>

namespace sbs {
namespace physics {
namespace xpbd {

// forward declare
class tetrahedral_mesh_t;

class green_constraint_t : public constraint_t
{
  public:
    using scalar_type       = typename constraint_t::scalar_type;
    using index_type        = std::uint32_t;
    using body_ptr_type     = tetrahedral_mesh_t*;
    using position_key_type = std::pair<body_ptr_type, index_type>;

    green_constraint_t(
        scalar_type const alpha,
        position_key_type const& vb1,
        position_key_type const& vb2,
        position_key_type const& vb3,
        position_key_type const& vb4,
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
    body_ptr_type b1_;
    index_type v1_;

    body_ptr_type b2_;
    index_type v2_;

    body_ptr_type b3_;
    index_type v3_;

    body_ptr_type b4_;
    index_type v4_;

    Eigen::Matrix3d DmInv_;
    scalar_type V0_;
    scalar_type mu_;
    scalar_type lambda_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H