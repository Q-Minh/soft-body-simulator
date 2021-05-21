#ifndef SBS_RENDERING_CAMERA_H
#define SBS_RENDERING_CAMERA_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace sbs {
namespace rendering {

class camera_t
{
  public:
    enum class movement_t { forward, backward, left, right, up, down };

    camera_t(
        glm::vec3 const& position = glm::vec3(0.f, 0.f, 5.f),
        glm::vec3 const& up       = glm::vec3(0.f, 1.f, 0.f),
        float yaw                 = -90.f,
        float pitch               = 0.f,
        float near                = 0.1f,
        float far                 = 100.f);

    glm::mat4 view_matrix() const { return glm::lookAt(position_, position_ + front_, up_); }
    glm::mat4 projection_matrix(float aspect_ratio) const
    {
        return glm::perspective(glm::radians(zoom_), aspect_ratio, near_, far_);
    }
    glm::vec3 const& position() const { return position_; }

    void handle_keyboard(movement_t direction, float dt);

    void handle_mouse_movement(float xoffset, float yoffset, bool constrain_pitch = true);

    void handle_mouse_scroll(float yoffset);

  private:
    void update_camera_vectors();

  private:
    glm::vec3 position_;
    glm::vec3 front_;
    glm::vec3 up_;
    glm::vec3 right_;
    glm::vec3 world_up_;

    float yaw_;
    float pitch_;
    float speed_;
    float mouse_sensitivity_;
    float zoom_;

    float near_;
    float far_;
};

} // namespace rendering
} // namespace sbs

#endif // SBS_RENDERING_CAMERA_H