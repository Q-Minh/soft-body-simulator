#include "sbs/rendering/trackball_rotation_adapter.h"

#include <Eigen/Geometry>

namespace sbs {
namespace rendering {
namespace detail {

static Eigen::Vector3d geometric_center(common::dynamic_surface_mesh const& mesh)
{
    Eigen::Vector3d mu{0., 0., 0.};

    for (std::size_t vi = 0u; vi < mesh.vertex_count(); ++vi)
    {
        auto const& v = mesh.vertex(vi);
        mu += v.position;
    }

    mu = mu.array() / mesh.vertex_count();
    return mu;
}

static void rotate(common::dynamic_surface_mesh& mesh, Eigen::Matrix3d const& rotation)
{
    for (std::size_t vi = 0u; vi < mesh.vertex_count(); ++vi)
    {
        auto& v            = mesh.mutable_vertex(vi);
        Eigen::Vector3d& p = v.position;
        p                  = rotation * p;
    }
}

static void translate(common::dynamic_surface_mesh& mesh, Eigen::Vector3d const& t)
{
    for (std::size_t vi = 0u; vi < mesh.vertex_count(); ++vi)
    {
        auto& v = mesh.mutable_vertex(vi);
        v.position += t;
    }
}

} // namespace detail

trackball_rotation_adapter_t::trackball_rotation_adapter_t(
    common::dynamic_surface_mesh* mesh,
    double rotation_speed,
    double translation_speed)
    : rotated_mesh_(mesh),
      rotation_speed_(rotation_speed),
      translation_speed_(translation_speed),
      pitch_angle_(0.),
      yaw_angle_(0.),
      pitch_axis_({1., 0., 0.}),
      yaw_axis_({0., 0., 1.})
{
}

double trackball_rotation_adapter_t::rotation_speed() const
{
    return rotation_speed_;
}

void trackball_rotation_adapter_t::set_rotation_speed(double speed)
{
    rotation_speed_ = speed;
}

double trackball_rotation_adapter_t::translation_speed() const
{
    return translation_speed_;
}

void trackball_rotation_adapter_t::set_translation_speed(double speed)
{
    translation_speed_ = speed;
}

void trackball_rotation_adapter_t::rotate(double dx, double dy)
{
    yaw_angle_   = rotation_speed_ * dx;
    pitch_angle_ = rotation_speed_ * dy;

    Eigen::AngleAxisd const yaw(yaw_angle_, yaw_axis_);
    Eigen::AngleAxisd const pitch(pitch_angle_, pitch_axis_);
    Eigen::Matrix3d const rotation = pitch.toRotationMatrix() * yaw.toRotationMatrix();

    Eigen::Vector3d const geometric_center = detail::geometric_center(*rotated_mesh_);
    detail::translate(*rotated_mesh_, -geometric_center);
    detail::rotate(*rotated_mesh_, rotation);
    detail::translate(*rotated_mesh_, geometric_center);
}

void trackball_rotation_adapter_t::translate(double dx, double dy)
{
    Eigen::Vector3d const t            = dx * pitch_axis_ + dy * yaw_axis_;
    Eigen::Vector3d const displacement = translation_speed_ * t;
    detail::translate(*rotated_mesh_, displacement);
}

void trackball_rotation_adapter_t::set_pitch_axis(Eigen::Vector3d const& pitch_axis)
{
    pitch_axis_ = pitch_axis.normalized();
}

void trackball_rotation_adapter_t::set_yaw_axis(Eigen::Vector3d const& yaw_axis)
{
    yaw_axis_ = yaw_axis.normalized();
}

common::dynamic_surface_mesh* trackball_rotation_adapter_t::mesh() const
{
    return rotated_mesh_;
}

common::dynamic_surface_mesh* trackball_rotation_adapter_t::mesh()
{
    return rotated_mesh_;
}

} // namespace rendering
} // namespace sbs