//#ifndef SBS_PHYSICS_CUTTING_VIRTUAL_SCALPEL_H
//#define SBS_PHYSICS_CUTTING_VIRTUAL_SCALPEL_H
//
//#include "common/mesh.h"
//#include "common/primitive.h"
//
//namespace sbs {
//namespace physics {
//namespace cutting {
//
//class virtual_scalpel_h
//{
//  public:
//    using self_type = virtual_scalpel_h;
//
//    virtual_scalpel_h(
//        common::line_segment_t const& l,
//        common::shared_vertex_surface_mesh_t const& render_model,
//        double rotation_speed = 1e-3,
//        double eps            = 1e-10);
//    virtual_scalpel_h(self_type const& rhs) = default;
//    virtual_scalpel_h(self_type&& rhs)      = default;
//
//    double rotation_speed() const;
//    void set_rotation_speed(double speed);
//    void rotate(double dx, double dy);
//    void set_translation(Eigen::Vector3d const& t);
//    void set_pitch_axis(Eigen::Vector3d const& pitch_axis);
//    void set_yaw_axis(Eigen::Vector3d const& yaw_axis);
//
//    void step();
//    common::shared_vertex_surface_mesh_t get_scalpel_render_model() const;
//    common::shared_vertex_surface_mesh_t const& get_render_swept_surface() const;
//    common::shared_vertex_surface_mesh_t const& get_swept_surface() const;
//    common::line_segment_t const& get_cutting_edge() const;
//    bool can_surface_be_swept() const;
//
//  protected:
//    enum class state_t { in_transition, sweeping };
//    void compute_swept_surface();
//
//  private:
//    common::line_segment_t l1_;
//    common::line_segment_t l2_;
//    common::shared_vertex_surface_mesh_t swept_surface_;
//    common::shared_vertex_surface_mesh_t last_valid_swept_surface_;
//    Eigen::Vector3d previous_translation_, translation_;
//
//    double rotation_speed_;
//    double previous_pitch_angle_, pitch_angle_, previous_yaw_angle_, yaw_angle_;
//    Eigen::Vector3d pitch_axis_, yaw_axis_;
//    double eps_;
//
//    common::shared_vertex_surface_mesh_t render_model_;
//    state_t state_;
//};
//
//} // namespace cutting
//} // namespace physics
//} // namespace sbs
//
//#endif // SBS_PHYSICS_CUTTING_VIRTUAL_SCALPEL_H