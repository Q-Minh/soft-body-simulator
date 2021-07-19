//#include "sbs/rendering/trackball_rotation_adapter.h"
//
//#include <Eigen/Geometry>
//
//namespace sbs {
//namespace rendering {
//
//trackball_rotation_adapter_t::trackball_rotation_adapter_t(
//    common::shared_vertex_surface_mesh_t const& mesh,
//    double rotation_speed)
//    : rotated_mesh_(mesh),
//      rotation_speed_(rotation_speed),
//      pitch_angle_(0.),
//      yaw_angle_(0.),
//      pitch_axis_({1., 0., 0.}),
//      yaw_axis_({0., 0., 1.})
//{
//}
//
//double trackball_rotation_adapter_t::rotation_speed() const
//{
//    return rotation_speed_;
//}
//
//void trackball_rotation_adapter_t::set_rotation_speed(double speed)
//{
//    rotation_speed_ = speed;
//}
//
//void trackball_rotation_adapter_t::rotate(double dx, double dy)
//{
//    yaw_angle_   = rotation_speed_ * dx;
//    pitch_angle_ = rotation_speed_ * dy;
//
//    Eigen::AngleAxisd const yaw(yaw_angle_, yaw_axis_);
//    Eigen::AngleAxisd const pitch(pitch_angle_, pitch_axis_);
//    Eigen::Matrix3d const rotation = pitch.toRotationMatrix() * yaw.toRotationMatrix();
//
//    Eigen::Vector3d const geometric_center = rotated_mesh_.vertices().rowwise().mean();
//    rotated_mesh_.vertices().colwise() -= geometric_center;
//    rotated_mesh_.vertices() = rotation * rotated_mesh_.vertices();
//    rotated_mesh_.vertices().colwise() += geometric_center;
//}
//
//void trackball_rotation_adapter_t::set_pitch_axis(Eigen::Vector3d const& pitch_axis)
//{
//    pitch_axis_ = pitch_axis.normalized();
//}
//
//void trackball_rotation_adapter_t::set_yaw_axis(Eigen::Vector3d const& yaw_axis)
//{
//    yaw_axis_ = yaw_axis.normalized();
//}
//
//common::shared_vertex_surface_mesh_t const& trackball_rotation_adapter_t::mesh() const
//{
//    return rotated_mesh_;
//}
//
//} // namespace rendering
//} // namespace sbs