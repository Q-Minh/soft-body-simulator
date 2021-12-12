#ifndef SBS_GEOMETRY_LINE_H
#define SBS_GEOMETRY_LINE_H

#include "sbs/aliases.h"
#include "sbs/common/mesh.h"

#include <Eigen/Core>
#include <array>
#include <functional>
#include <optional>
#include <utility>
#include <vector>

namespace sbs {
namespace geometry {

template <unsigned int Dimension>
class line_segment_t
{
  public:
    using point_type  = Eigen::Vector<scalar_type, Dimension>;
    using vector_type = Eigen::Vector<scalar_type, Dimension>;
    using self_type   = line_segment_t<Dimension>;

    line_segment_t(point_type const& p1, point_type const& p2) : p_(p1), v_(), length_()
    {
        vector_type const p12 = p2 - p1;
        length_               = p12.norm();
        v_                    = p12 / length_;
    }
    line_segment_t(point_type const& p, vector_type const& v, scalar_type length)
        : p_(p), v_(v), length_(length)
    {
    }

    point_type p1() const { return p_; }
    point_type p2() const { return p_ + length_ * v_; }

    vector_type direction() const { return v_; }
    scalar_type length() const { return length_; }

    /**
     * @brief Triangulated swept surface from this line segment to the other line segment.
     * Discretizes the swept surface with 2 triangles only. Orientation of the surface goes
     * in the direction of p1() to p2().
     * @return
     */
    std::pair<std::array<point_type, 4>, std::array<index_type, 6>>
    swept_surface(self_type const& other) const;

  private:
    point_type p_;       ///< starting point on the line segment
    vector_type v_;      ///< direction of the line
    scalar_type length_; ///< length of the line segment
};

template <unsigned int Dimension>
inline std::
    pair<std::array<typename line_segment_t<Dimension>::point_type, 4>, std::array<index_type, 6>>
    line_segment_t<Dimension>::swept_surface(self_type const& other) const
{
    std::array<point_type, 4> points{};
    std::array<index_type, 6> indices{};

    points[0] = this->p1();
    points[1] = this->p2();
    points[2] = other.p1();
    points[3] = other.p2();

    indices[0] = 0u;
    indices[1] = 1u;
    indices[2] = 2u;

    indices[3] = 2u;
    indices[4] = 1u;
    indices[5] = 3u;

    return std::make_pair(points, indices);
}

std::tuple<std::vector<Eigen::Vector3d>, std::vector<index_type>, line_segment_t<3>> swept_surface(
    line_segment_t<3> const& line,
    std::function<std::optional<Eigen::Vector3d>(scalar_type)> const& parametric_curve,
    scalar_type t0,
    scalar_type dt,
    index_type n);

class swept_line_segment_t : public common::shared_vertex_surface_mesh_i
{
  public:
    swept_line_segment_t(line_segment_t<3> const& starting_line_segment);

    void
    set_parametric_curve(std::function<std::optional<Eigen::Vector3d>(scalar_type)> const& curve)
    {
        parametric_curve_ = curve;
    }
    void set_dt(scalar_type dt) { dt_ = dt; }
    void sweep(scalar_type delta_t);

    void set_color(Eigen::Vector3f const& c) { color_ = c; }

    virtual std::size_t triangle_count() const { return triangles_.size(); };
    virtual std::size_t vertex_count() const { return vertices_.size(); };

    virtual vertex_type vertex(std::size_t vi) const { return vertices_[vi]; };
    virtual triangle_type triangle(std::size_t f) const { return triangles_[f]; };

    virtual void prepare_vertices_for_rendering() override;
    virtual void prepare_indices_for_rendering() override;

  private:
    line_segment_t<3> starting_line_segment_;
    scalar_type t_;
    scalar_type dt_;
    std::function<std::optional<Eigen::Vector3d>(scalar_type)> parametric_curve_;

    std::vector<vertex_type> vertices_;
    std::vector<triangle_type> triangles_;

    Eigen::Vector3f color_;
};

} // namespace geometry
} // namespace sbs

#endif // SBS_GEOMETRY_LINE_H
