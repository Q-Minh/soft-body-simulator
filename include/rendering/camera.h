#ifndef SBS_RENDERING_CAMERA_H
#define SBS_RENDERING_CAMERA_H

#include <Eigen/Core>
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

    glm::mat4 view_gl() const { return glm::lookAt(position_, position_ + front_, up_); }
    glm::mat4 projection_gl(float aspect_ratio) const
    {
        return glm::perspective(glm::radians(zoom_), aspect_ratio, near_, far_);
    }
    Eigen::Matrix4d view_matrix() const
    {
        glm::mat4 const view = view_gl();
        Eigen::Matrix4d view_mat{};
        view_mat.col(0) = Eigen::Vector4d{view[0].x, view[0].y, view[0].z, view[0].w};
        view_mat.col(1) = Eigen::Vector4d{view[1].x, view[1].y, view[1].z, view[1].w};
        view_mat.col(2) = Eigen::Vector4d{view[2].x, view[2].y, view[2].z, view[2].w};
        view_mat.col(3) = Eigen::Vector4d{view[3].x, view[3].y, view[3].z, view[3].w};
        return view_mat;
    }
    Eigen::Matrix4d projection_matrix(float aspect_ratio) const
    {
        glm::mat4 const proj = projection_gl(aspect_ratio);
        Eigen::Matrix4d proj_mat{};
        proj_mat.col(0) = Eigen::Vector4d{proj[0].x, proj[0].y, proj[0].z, proj[0].w};
        proj_mat.col(1) = Eigen::Vector4d{proj[1].x, proj[1].y, proj[1].z, proj[1].w};
        proj_mat.col(2) = Eigen::Vector4d{proj[2].x, proj[2].y, proj[2].z, proj[2].w};
        proj_mat.col(3) = Eigen::Vector4d{proj[3].x, proj[3].y, proj[3].z, proj[3].w};
        return proj_mat;
    }

    glm::vec3 front() const { return front_; }
    glm::vec3 up() const { return up_; }
    glm::vec3 right() const { return right_; }

    glm::vec3 const& position() const { return position_; }

    void handle_keyboard(movement_t direction, float dt);

    void handle_mouse_movement(float xoffset, float yoffset, bool constrain_pitch = true);

    void handle_mouse_scroll(float yoffset);

    void set_mouse_sensitivity(float sensitivity);

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