#include "sbs/rendering/physics_timestep_throttler.h"

#include <chrono>

namespace sbs {
namespace rendering {

physics_timestep_throttler_t::physics_timestep_throttler_t(
    double timestep,
    std::function<void(
        physics_timestep_throttler_t& /*throttler*/,
        double /*physics_dt*/,
        common::scene_t& /*scene*/)> step)
    : tb_{0.}, timestep_(timestep), fps_{0u}, step_(step), are_physics_active_{false}
{
}

void physics_timestep_throttler_t::operator()(double frame_dt, common::scene_t& scene)
{
    tb_ += frame_dt;
    auto const time_between_frames = tb_;

    // If the elapsed time between the last physics update was less than
    // the physics timestep, we don't update physics.
    if (tb_ < timestep_)
        return;

    double const num_timesteps_elapsed = std::floor(tb_ / timestep_);
    tb_ -= num_timesteps_elapsed * timestep_;

    if (!are_physics_active_)
        return;

    auto const begin = std::chrono::steady_clock::now();

    step_(*this, timestep_, scene);

    auto const end      = std::chrono::steady_clock::now();
    auto const duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();

    fps_ = static_cast<std::size_t>(1e9 / static_cast<double>(duration));
}

std::size_t physics_timestep_throttler_t::fps() const
{
    return fps_;
}

double physics_timestep_throttler_t::timestep() const
{
    return timestep_;
}
double& physics_timestep_throttler_t::timestep()
{
    return timestep_;
}

bool physics_timestep_throttler_t::are_physics_active() const
{
    return are_physics_active_;
}

void physics_timestep_throttler_t::deactivate_physics()
{
    are_physics_active_ = false;
}

void physics_timestep_throttler_t::activate_physics()
{
    are_physics_active_ = true;
}

} // namespace rendering
} // namespace sbs