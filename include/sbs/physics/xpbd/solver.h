#ifndef SBS_PHYSICS_XPBD_SOLVER_H
#define SBS_PHYSICS_XPBD_SOLVER_H

#include "sbs/aliases.h"

#include <cstddef>

namespace sbs {
namespace physics {
namespace xpbd {

// Forward declares
class simulation_t;

class solver_t
{
  public:
    virtual void solve(simulation_t& simulation, scalar_type dt, std::size_t iterations) = 0;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_SOLVER_H