#ifndef SBS_PHYSICS_XPBD_DISTANCE_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_DISTANCE_CONSTRAINT_H

#include "constraint.h"

namespace sbs {
namespace physics {
namespace xpbd {

class distance_constraint_t : public constraint_t
{
  public:
    using scalar_type     = typename constraint_t::scalar_type;
    using index_type      = typename constraint_t::index_type;
    using index_pair_type = std::pair<index_type, index_type>;
    using positions_type  = typename constraint_t::positions_type;
    using masses_type     = typename constraint_t::masses_type;

    distance_constraint_t(
        scalar_type const alpha,
        positions_type const& positions,
        index_pair_type const& vb1,
        index_pair_type const& vb2);

    virtual void project(
        std::vector<positions_type>& positions,
        std::vector<masses_type> const& masses,
        scalar_type& lagrange_multiplier,
        scalar_type const dt) const override;

    scalar_type evaluate(Eigen::Vector3d const& p1, Eigen::Vector3d const& p2) const;

  private:
    index_type v1_;
    index_type v2_;
    index_type b1_;
    index_type b2_;
    scalar_type d_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_DISTANCE_CONSTRAINT_H