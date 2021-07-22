#ifndef SBS_RENDERING_TRACKBALL_ROTATION_ADAPTER_H
#define SBS_RENDERING_TRACKBALL_ROTATION_ADAPTER_H

#include "sbs/common/mesh.h"

#include <Eigen/Core>

namespace sbs {
namespace rendering {

class trackball_rotation_adapter_t
{
  public:
    trackball_rotation_adapter_t() = default;
    trackball_rotation_adapter_t(
        common::dynamic_surface_mesh* mesh,
        double rotation_speed    = 1e-3,
        double translation_speed = 1e-1);

    double rotation_speed() const;
    void set_rotation_speed(double speed);
    double translation_speed() const;
    void set_translation_speed(double speed);

    void rotate(double dx, double dy);
    void translate(double dx, double dy);
    void set_pitch_axis(Eigen::Vector3d const& pitch_axis);
    void set_yaw_axis(Eigen::Vector3d const& yaw_axis);

    common::dynamic_surface_mesh* mesh() const;
    common::dynamic_surface_mesh* mesh();

  private:
    common::dynamic_surface_mesh* rotated_mesh_;

    double rotation_speed_;
    double translation_speed_;
    double pitch_angle_, yaw_angle_;
    Eigen::Vector3d pitch_axis_, yaw_axis_;
};

} // namespace rendering
} // namespace sbs

#endif // SBS_RENDERING_TRACKBALL_ROTATION_ADAPTER_H