#ifndef SBS_PHYSICS_XPBD_SOLVER_H
#define SBS_PHYSICS_XPBD_SOLVER_H

#include "constraint.h"

#include <memory>
#include <vector>

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
    double young_modulus = 1e9;
    double poisson_ratio = 0.48;
};

class solver_t
{
  public:
    void setup(
        std::vector<std::shared_ptr<common::node_t>>* bodies,
        std::vector<simulation_parameters_t> const& per_body_simulation_parameters);

    void step(double timestep, std::uint32_t iterations, std::uint32_t substeps);
    void notify_topology_changed();

    simulation_parameters_t const& get_simulation_parameters_for_body(std::uint32_t index) const;
    simulation_parameters_t& get_simulation_parameters_for_body(std::uint32_t index);

  protected:
    void handle_collisions(std::vector<Eigen::Matrix3Xd> const& P);

  private:
    std::vector<std::shared_ptr<common::node_t>>* bodies_;
    std::vector<simulation_parameters_t> per_body_simulation_parameters_;

    std::vector<std::unique_ptr<constraint_t>> constraints_;
    std::vector<std::unique_ptr<constraint_t>> collision_constraints_;
    // std::vector<std::uint32_t>
    //    garbage_collected_constraints_; ///< We collect, in this list, indices of constraints that
    //                                    ///< have been removed, because it is too expensive to
    //                                    ///< remove them directly from the list of constraints
    //                                    ///< or because we wish to preserve the ordering of the
    //                                    ///< constraints
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_XPBD_SOLVER_H