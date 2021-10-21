#ifndef SBS_PHYSICS_XPBD_SIMULATION_PARAMETERS_H
#define SBS_PHYSICS_XPBD_SIMULATION_PARAMETERS_H

#include <sbs/aliases.h>

namespace sbs {
namespace physics {
namespace xpbd {

struct simulation_parameters_t
{
    scalar_type compliance = 1e-4;
    scalar_type damping    = 0.;

    /**
     * FEM parameters
     */
    scalar_type young_modulus = 1e6;
    scalar_type poisson_ratio = 0.30;

    /*
     * Spring parameters
     */
    scalar_type hooke_coefficient = 1.;

    scalar_type collision_compliance = 1e-8;
    scalar_type collision_damping    = 0.;

    /**
     * Meshless parameters
     */
    scalar_type positional_penalty_strength = 1.;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_SIMULATION_PARAMETERS_H