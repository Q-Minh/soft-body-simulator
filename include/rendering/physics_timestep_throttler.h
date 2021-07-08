#ifndef SBS_RENDERING_PHYSICS_TIMESTEP_THROTTLER_H
#define SBS_RENDERING_PHYSICS_TIMESTEP_THROTTLER_H

#include <functional>

namespace sbs {

// forward declares
namespace common {

struct scene_t;

} // namespace common

namespace rendering {

class physics_timestep_throttler_t
{
  public:
    physics_timestep_throttler_t(
        double timestep,
        std::function<void(
            physics_timestep_throttler_t& /*throttler*/,
            double /*physics_dt*/,
            sbs::common::scene_t& /*scene*/)> step);
    physics_timestep_throttler_t(physics_timestep_throttler_t const&) = default;
    physics_timestep_throttler_t(physics_timestep_throttler_t&&)      = default;

    void operator()(double frame_dt, sbs::common::scene_t& scene);

    std::size_t fps() const;

    double timestep() const;
    double& timestep();

  private:
    double tb_;       ///< Time elapsed since last physics step
    double timestep_; ///< Desired timestep
    std::size_t fps_; ///< Computed frames per second
    std::function<
        void(physics_timestep_throttler_t&, double /*physics_dt*/, sbs::common::scene_t& /*scene*/)>
        step_; ///< Physics step to execute with throttling
};

} // namespace rendering
} // namespace sbs

#endif // SBS_RENDERING_PHYSICS_TIMESTEP_THROTTLER_H