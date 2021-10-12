#ifndef SBS_PHYSICS_TIMESTEP_H
#define SBS_PHYSICS_TIMESTEP_H

#include <cstddef>
#include <memory>
#include <sbs/aliases.h>

namespace sbs {
namespace physics {

class simulation_t;
class solver_t;

class timestep_t
{
  public:
    void step(simulation_t& simulation);

    scalar_type dt() const;
    scalar_type& dt();

    std::size_t iterations() const;
    std::size_t& iterations();

    std::size_t substeps() const;
    std::size_t& substeps();

    std::unique_ptr<solver_t> const& solver() const;
    std::unique_ptr<solver_t>& solver();

  private:
    scalar_type dt_;
    std::size_t iterations_;
    std::size_t substeps_;
    std::unique_ptr<solver_t> solver_;
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_TIMESTEP_H