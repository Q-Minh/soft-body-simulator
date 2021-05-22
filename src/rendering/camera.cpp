#include "rendering/camera.h"

namespace sbs {
namespace rendering {

camera_t::camera_t(
    glm::vec3 const& position,
    glm::vec3 const& up,
    float yaw,
    float pitch,
    float near,
    float far)
    : position_(position),
      front_(glm::vec3(0.f, 0.f, -1.f)),
      up_(),
      right_(),
      world_up_(up),
      yaw_(yaw),
      pitch_(pitch),
      speed_(2.5f),
      mouse_sensitivity_(0.1f),
      zoom_(45.f),
      near_(near),
      far_(far)
{
    update_camera_vectors();
}

void camera_t::handle_keyboard(movement_t direction, float dt)
{
    float velocity = speed_ * dt;
    if (direction == movement_t::forward)
        position_ += front_ * velocity;
    if (direction == movement_t::backward)
        position_ -= front_ * velocity;
    if (direction == movement_t::left)
        position_ -= right_ * velocity;
    if (direction == movement_t::right)
        position_ += right_ * velocity;
    if (direction == movement_t::up)
        position_ += up_ * velocity;
    if (direction == movement_t::down)
        position_ -= up_ * velocity;
}

void camera_t::handle_mouse_movement(float xoffset, float yoffset, bool constrain_pitch)
{
    xoffset *= mouse_sensitivity_;
    yoffset *= mouse_sensitivity_;

    yaw_ += xoffset;
    pitch_ += yoffset;

    if (constrain_pitch)
    {
        if (pitch_ > 89.0f)
            pitch_ = 89.0f;
        if (pitch_ < -89.0f)
            pitch_ = -89.0f;
    }

    update_camera_vectors();
}

void camera_t::handle_mouse_scroll(float yoffset)
{
    zoom_ -= (float)yoffset;
    if (zoom_ < 1.0f)
        zoom_ = 1.0f;
    if (zoom_ > 45.0f)
        zoom_ = 45.0f;
}

void camera_t::update_camera_vectors()
{
    glm::vec3 front;
    front.x = cos(glm::radians(yaw_)) * cos(glm::radians(pitch_));
    front.y = sin(glm::radians(pitch_));
    front.z = sin(glm::radians(yaw_)) * cos(glm::radians(pitch_));
    front_  = glm::normalize(front);
    right_  = glm::normalize(glm::cross(front_, world_up_));
    up_     = glm::normalize(glm::cross(right_, front_));
}

void camera_t::set_mouse_sensitivity(float sensitivity)
{
    mouse_sensitivity_ = sensitivity;
}

} // namespace rendering
} // namespace sbs