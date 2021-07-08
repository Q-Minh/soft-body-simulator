#ifndef SBS_PHYSICS_XPBD_SIMULATION_PARAMETERS_H
#define SBS_PHYSICS_XPBD_SIMULATION_PARAMETERS_H

namespace sbs {
namespace physics {
namespace xpbd {

enum class constraint_type_t { green, distance };

struct simulation_parameters_t
{
    double alpha                      = 0.1;
    constraint_type_t constraint_type = constraint_type_t::green;

    /**
     * FEM parameters
     */
    double young_modulus = 1e6;
    double poisson_ratio = 0.45;

    /*
     * Spring parameters
     */
    double hooke_coefficient = 1.;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_SIMULATION_PARAMETERS_H