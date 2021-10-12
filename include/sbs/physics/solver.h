#ifndef SBS_PHYSICS_SOLVER_H
#define SBS_PHYSICS_SOLVER_H

#include <cstddef>
#include <sbs/aliases.h>

namespace sbs {
namespace physics {

// Forward declares
class simulation_t;

class solver_t
{
  public:
    virtual void solve(simulation_t& simulation, scalar_type dt, std::size_t iterations) = 0;
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_SOLVER_H