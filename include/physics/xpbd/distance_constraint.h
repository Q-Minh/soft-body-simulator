#ifndef SBS_PHYSICS_XPBD_DISTANCE_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_DISTANCE_CONSTRAINT_H

#include "constraint.h"

#include <Eigen/Core>

namespace sbs {
namespace physics {
namespace xpbd {

class distance_constraint_t : public constraint_t
{
  public:
    using scalar_type       = typename constraint_t::scalar_type;
    using index_type        = std::uint32_t;
    using body_ptr_type     = tetrahedral_mesh_t*;
    using position_key_type = std::pair<body_ptr_type, index_type>;

    distance_constraint_t(
        scalar_type const alpha,
        position_key_type const& vb1,
        position_key_type const& vb2);

    virtual void project(
        std::vector<std::shared_ptr<xpbd::tetrahedral_mesh_t>> const& bodies,
        scalar_type& lagrange_multiplier,
        scalar_type const dt) const override;

    scalar_type evaluate(Eigen::Vector3d const& p1, Eigen::Vector3d const& p2) const;

  private:
    body_ptr_type b1_;
    index_type v1_;
    body_ptr_type b2_;
    index_type v2_;
    scalar_type d_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_DISTANCE_CONSTRAINT_H