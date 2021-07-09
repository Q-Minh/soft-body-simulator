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

    /**
     * @brief Handler for the renderer's loop's new_physics_update event
     * @param frame_dt Renderer's duration/elapsed time since last frame
     * @param scene The renderer's scene
     */
    void operator()(double frame_dt, sbs::common::scene_t& scene);

    /**
     * @brief Gets the estimated fps of the physics simulation loop. Does not correspond to the
     * rendering loop's fps. Also does not correspond to inverse of the timestep. This fps
     * measure is really computing the time it takes for physics steps to complete.
     * @return Estimated physics fps.
     */
    std::size_t fps() const;

    /**
     * @brief The desired timestep of the physics loop. The throttler will
     * do its best to only run the physics loop at the desired time step.
     * @return The current desired timestep of the physics loop.
     */
    double timestep() const;
    double& timestep();

    bool are_physics_active() const;
    void deactivate_physics();
    void activate_physics();

  private:
    double tb_;       ///< Time elapsed since last physics step
    double timestep_; ///< Desired timestep
    std::size_t fps_; ///< Computed frames per second
    std::function<
        void(physics_timestep_throttler_t&, double /*physics_dt*/, sbs::common::scene_t& /*scene*/)>
        step_; ///< Physics step to execute with throttling
    bool are_physics_active_;
};

} // namespace rendering
} // namespace sbs

#endif // SBS_RENDERING_PHYSICS_TIMESTEP_THROTTLER_H