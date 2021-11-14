#ifndef SBS_PHYSICS_XPBD_CONSTRAINT_H
#define SBS_PHYSICS_XPBD_CONSTRAINT_H

#include "sbs/aliases.h"

namespace sbs {
namespace physics {
namespace xpbd {

// Forward declares
class simulation_t;

class constraint_t
{
  public:
    constraint_t(scalar_type alpha, scalar_type beta);

    void prepare_for_projection(simulation_t& simulation);
    virtual void project_positions(simulation_t& simulation, scalar_type dt) = 0;

    virtual ~constraint_t() = default;

    scalar_type alpha() const;
    scalar_type beta() const;
    scalar_type lambda() const;

    scalar_type compliance() const;
    scalar_type damping() const;

  protected:
    virtual void prepare_for_projection_impl(simulation_t& simulation) {}

    scalar_type alpha_;
    scalar_type beta_;
    scalar_type lagrange_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_CONSTRAINT_H