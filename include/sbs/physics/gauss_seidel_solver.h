#ifndef SBS_PHYSICS_GAUSS_SEIDEL_SOLVER_H
#define SBS_PHYSICS_GAUSS_SEIDEL_SOLVER_H

#include <sbs/physics/solver.h>

namespace sbs {
namespace physics {

class gauss_seidel_solver_t : public solver_t
{
  public:
    virtual void solve(simulation_t& simulation, scalar_type dt, std::size_t iterations) override;
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_GAUSS_SEIDEL_SOLVER_H