#ifndef SBS_PHYSICS_SIMULATION_H
#define SBS_PHYSICS_SIMULATION_H

#include "sbs/aliases.h"
#include "sbs/physics/body/body.h"
#include "sbs/physics/collision/cd_system.h"
#include "sbs/physics/xpbd/particle.h"
#include "sbs/physics/xpbd/constraint.h"
#include "sbs/physics/xpbd/simulation_parameters.h"

#include <vector>

namespace sbs {
namespace physics {
namespace xpbd {

class simulation_t
{
  public:
    void use_collision_detection_system(std::unique_ptr<collision::cd_system_t> cd_system);
    void add_particle(particle_t const& p, index_type const body_idx);
    index_type add_body(std::unique_ptr<body::body_t> body);
    index_type add_body();
    void add_constraint(std::unique_ptr<xpbd::constraint_t> constraint);
    void remove_constraint(index_type const constraint_idx);
    void add_collision_constraint(std::unique_ptr<xpbd::constraint_t> collision_constraint);

    std::vector<std::vector<particle_t>> const& particles() const;
    std::vector<std::vector<particle_t>>& particles();
    std::vector<std::unique_ptr<body::body_t>> const& bodies() const;
    std::vector<std::unique_ptr<body::body_t>>& bodies();
    std::vector<std::unique_ptr<xpbd::constraint_t>> const& constraints() const;
    std::vector<std::unique_ptr<xpbd::constraint_t>>& constraints();
    std::vector<std::unique_ptr<xpbd::constraint_t>> const& collision_constraints() const;
    std::vector<std::unique_ptr<xpbd::constraint_t>>& collision_constraints();
    std::unique_ptr<collision::cd_system_t> const& collision_detection_system() const;
    xpbd::simulation_parameters_t const& simulation_parameters() const;
    xpbd::simulation_parameters_t& simulation_parameters();

  private:
    std::vector<std::vector<particle_t>> particles_;
    std::vector<std::unique_ptr<body::body_t>> bodies_;
    std::vector<std::unique_ptr<xpbd::constraint_t>> constraints_;
    std::vector<std::unique_ptr<xpbd::constraint_t>> collision_constraints_;
    std::unique_ptr<collision::cd_system_t> cd_system_;
    xpbd::simulation_parameters_t simulation_parameters_;
};

} // namespace xpbd
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_SIMULATION_H