#include "physics\cutting\virtual_scalpel.h"

#include <Eigen/Geometry>

namespace sbs {
namespace physics {
namespace cutting {

virtual_scalpel_h::virtual_scalpel_h(
    common::line_segment_t const& l,
    common::shared_vertex_surface_mesh_t const& render_model,
    double rotation_speed,
    double eps)
    : l1_(l),
      l2_(l),
      swept_surface_(),
      previous_translation_{0., 0., 0.},
      translation_{0., 0., 0.},
      rotation_speed_(rotation_speed),
      previous_pitch_angle_(0.),
      pitch_angle_(0.),
      previous_yaw_angle_(0.),
      yaw_angle_(0.),
      pitch_axis_(),
      yaw_axis_(),
      eps_(eps),
      render_model_(render_model),
      last_valid_swept_surface_(),
      state_(state_t::in_transition)
{
}

double virtual_scalpel_h::rotation_speed() const
{
    return rotation_speed_;
}

void virtual_scalpel_h::set_rotation_speed(double speed)
{
    rotation_speed_ = speed;
}

void virtual_scalpel_h::rotate(double dx, double dy)
{
    yaw_angle_   = rotation_speed_ * dx;
    pitch_angle_ = rotation_speed_ * dy;

    Eigen::AngleAxisd const yaw(yaw_angle_, yaw_axis_);
    Eigen::AngleAxisd const pitch(pitch_angle_, pitch_axis_);
    Eigen::Matrix3d const rotation = pitch.toRotationMatrix() * yaw.toRotationMatrix();

    // update render model
    Eigen::Vector3d geometric_mean = /*render_model_.vertices().rowwise().mean()*/ l2_.p;
    render_model_.vertices().colwise() -= geometric_mean;
    render_model_.vertices() = rotation * render_model_.vertices();
    render_model_.vertices().colwise() += geometric_mean;

    // update physical model
    geometric_mean = l2_.p;
    l2_.q          = rotation * (l2_.q - geometric_mean) + geometric_mean;

    double constexpr minimum_edge_length = 1e-10;
    bool const surface_can_be_swept      = (l2_.p - l1_.p).norm() > minimum_edge_length ||
                                      (l2_.q - l1_.q).norm() > minimum_edge_length;

    if (surface_can_be_swept && state_ == state_t::in_transition)
    {
        state_ = state_t::sweeping;
    }
}

void virtual_scalpel_h::set_translation(Eigen::Vector3d const& t)
{
    Eigen::Vector3d const displacement = t - translation_;
    translation_                       = t;

    // update render model
    render_model_.vertices().colwise() += displacement;
    // update physical model
    l2_.p += displacement;
    l2_.q += displacement;

    if (can_surface_be_swept() && state_ == state_t::in_transition)
    {
        state_ = state_t::sweeping;
    }
}

void virtual_scalpel_h::set_pitch_axis(Eigen::Vector3d const& pitch_axis)
{
    pitch_axis_ = pitch_axis;
}

void virtual_scalpel_h::set_yaw_axis(Eigen::Vector3d const& yaw_axis)
{
    yaw_axis_ = yaw_axis;
}

void virtual_scalpel_h::step()
{
    if (state_ == state_t::in_transition)
        return;

    compute_swept_surface();

    l1_ = l2_;

    state_ = state_t::in_transition;
}

common::shared_vertex_surface_mesh_t virtual_scalpel_h::get_scalpel_render_model() const
{
    common::shared_vertex_surface_mesh_t render_model{render_model_};
    render_model.set_color({0.f, 0.f, 0.f});
    Eigen::Vector3f const red{1.f, 0.f, 0.f};
    Eigen::Vector3f const blue{0.f, 0.f, 1.f};
    render_model.colors().col(0) = red;
    render_model.colors().col(1) = red;
    render_model.colors().col(2) = blue;
    render_model.colors().col(3) = blue;

    return render_model;
}

common::shared_vertex_surface_mesh_t const& virtual_scalpel_h::get_render_swept_surface() const
{
    return last_valid_swept_surface_;
}

common::shared_vertex_surface_mesh_t const& virtual_scalpel_h::get_swept_surface() const
{
    return swept_surface_;
}

common::line_segment_t const& virtual_scalpel_h::get_cutting_edge() const
{
    return l2_;
}

bool virtual_scalpel_h::can_surface_be_swept() const
{
    double constexpr minimum_edge_length = 1e-10;
    bool const surface_can_be_swept      = (l2_.p - l1_.p).norm() > minimum_edge_length ||
                                      (l2_.q - l1_.q).norm() > minimum_edge_length;
    return surface_can_be_swept;
}

void virtual_scalpel_h::compute_swept_surface()
{
    auto const& p1 = l1_.p;
    auto const& p2 = l1_.q;
    auto const& p3 = l2_.p;
    auto const& p4 = l2_.q;

    swept_surface_.vertices().resize(3, 4);
    swept_surface_.normals().resize(3, 4);
    swept_surface_.triangles().resize(3, 2);

    swept_surface_.vertices().col(0) = p1;
    swept_surface_.vertices().col(1) = p2;
    swept_surface_.vertices().col(2) = p3;
    swept_surface_.vertices().col(3) = p4;

    swept_surface_.triangles().col(0) =
        common::shared_vertex_surface_mesh_t::triangle_type{0, 1, 2};
    swept_surface_.triangles().col(1) =
        common::shared_vertex_surface_mesh_t::triangle_type{2, 1, 3};

    common::normal_t const n1 = (p2 - p1).cross(p3 - p1).normalized();
    common::normal_t const n2 = (p1 - p2).cross(p3 - p2).normalized();

    swept_surface_.normals().col(0) = n1;
    swept_surface_.normals().col(1) = ((n1 + n2) / 2.).normalized();
    swept_surface_.normals().col(2) = ((n1 + n2) / 2.).normalized();
    swept_surface_.normals().col(3) = n2;

    Eigen::Vector3f const red{1.f, 0.f, 0.f};
    Eigen::Vector3f const blue{0.f, 0.f, 1.f};
    swept_surface_.set_color(red);
    swept_surface_.colors().col(2) = blue;
    swept_surface_.colors().col(3) = blue;

    last_valid_swept_surface_ = swept_surface_;
}

} // namespace cutting
} // namespace physics
} // namespace sbs
