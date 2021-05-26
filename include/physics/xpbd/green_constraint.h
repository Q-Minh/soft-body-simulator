#ifndef SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H

#include "constraint.h"

namespace sbs {
namespace physics {
namespace xpbd {

class green_constraint_t : public constraint_t
{
  public:
    using scalar_type     = typename constraint_t::scalar_type;
    using index_type      = typename constraint_t::index_type;
    using index_pair_type = std::pair<index_type, index_type>;
    using positions_type  = typename constraint_t::positions_type;
    using masses_type     = typename constraint_t::masses_type;

    green_constraint_t(
        scalar_type const alpha,
        positions_type const& positions,
        index_pair_type const& vb1,
        index_pair_type const& vb2,
        index_pair_type const& vb3,
        index_pair_type const& vb4,
        scalar_type young_modulus,
        scalar_type poisson_ratio);

    virtual void project(
        std::vector<positions_type>& positions,
        std::vector<masses_type> const& masses,
        scalar_type& lagrange_multiplier,
        scalar_type const dt) const override;

  protected:
    scalar_type signed_volume(
        Eigen::Vector3d const& p1,
        Eigen::Vector3d const& p2,
        Eigen::Vector3d const& p3,
        Eigen::Vector3d const& p4) const;

  private:
    index_type v1_;
    index_type v2_;
    index_type v3_;
    index_type v4_;

    index_type b1_;
    index_type b2_;
    index_type b3_;
    index_type b4_;

    Eigen::Matrix3d DmInv_;
    scalar_type V0_;
    scalar_type mu_;
    scalar_type lambda_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H