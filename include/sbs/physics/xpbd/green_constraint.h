#ifndef SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_GREEN_CONSTRAINT_H

#include <Eigen/Core>
#include <sbs/physics/constraint.h>

namespace sbs {
namespace physics {

// Forward declares
class simulation_t;

namespace xpbd {

class green_constraint_t : public constraint_t
{
  public:
    green_constraint_t(
        scalar_type const alpha,
        scalar_type const beta,
        simulation_t const& simulation,
        index_type bi,
        index_type v1,
        index_type v2,
        index_type v3,
        index_type v4,
        scalar_type young_modulus,
        scalar_type poisson_ratio);

    virtual void project_positions(simulation_t& simulation, scalar_type dt) override;

  protected:
    scalar_type signed_volume(
        Eigen::Vector3d const& p1,
        Eigen::Vector3d const& p2,
        Eigen::Vector3d const& p3,
        Eigen::Vector3d const& p4) const;

  private:
    index_type bi_;
    index_type v1_;
    index_type v2_;
    index_type v3_;
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