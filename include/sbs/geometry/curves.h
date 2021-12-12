#ifndef SBS_GEOMETRY_CURVES_H
#define SBS_GEOMETRY_CURVES_H

#include "sbs/aliases.h"

#include <Eigen/Core>
#include <array>
#include <optional>
#include <vector>

namespace sbs {
namespace geometry {

template <unsigned int Dimension, unsigned int Order>
class bezier_curve_t;

/**
 * @brief Linear bezier curve
 */
template <unsigned int Dimension>
class bezier_curve_t<Dimension, 1>
{
  public:
    using point_type = Eigen::Vector<scalar_type, Dimension>;

    Eigen::Vector3d eval(scalar_type t) const
    {
        assert(t >= 0. && t <= 1.);
        return (1. - t) * P_[0] + t * P_[1];
    }
    std::array<point_type, 4> const& get_control_points() const { return P_; }
    std::array<point_type, 4>& get_control_points() { return P_; }

    void set_control_points(std::array<point_type, 2> const& P) { P_ = P; }

  private:
    std::array<point_type, 2> P_;
};

/**
 * @brief Quadratic bezier curve
 */
template <unsigned int Dimension>
class bezier_curve_t<Dimension, 2>
{
  public:
    using point_type = Eigen::Vector<scalar_type, Dimension>;

    Eigen::Vector3d eval(scalar_type t) const
    {
        assert(t >= 0. && t <= 1.);
        scalar_type const one_minus_t = 1. - t;
        scalar_type const c1          = one_minus_t * one_minus_t;
        scalar_type const c2          = 2. * one_minus_t * t;
        scalar_type const c3          = t * t;
        return c1 * P_[0] + c2 * P_[1] + c3 * P_[2];
    }
    std::array<point_type, 4> const& get_control_points() const { return P_; }
    std::array<point_type, 4>& get_control_points() { return P_; }

    void set_control_points(std::array<point_type, 3> const& P) { P_ = P; }

  private:
    std::array<point_type, 3> P_;
};

/**
 * @brief Cubic bezier curve
 */
template <unsigned int Dimension>
class bezier_curve_t<Dimension, 3>
{
  public:
    using point_type = Eigen::Vector<scalar_type, Dimension>;

    Eigen::Vector3d eval(scalar_type t) const
    {
        assert(t >= 0. && t <= 1.);
        scalar_type const one_minus_t   = 1. - t;
        scalar_type const one_minus_t_2 = one_minus_t * one_minus_t;
        scalar_type const one_minus_t_3 = one_minus_t_2 * one_minus_t;
        scalar_type const t2            = t * t;
        scalar_type const t3            = t2 * t;
        scalar_type const c1            = one_minus_t_3;
        scalar_type const c2            = 3. * one_minus_t_2 * t;
        scalar_type const c3            = 3. * one_minus_t * t2;
        scalar_type const c4            = t3;
        return c1 * P_[0] + c2 * P_[1] + c3 * P_[2] + c4 * P_[3];
    }
    std::array<point_type, 4> const& get_control_points() const { return P_; }
    std::array<point_type, 4>& get_control_points() { return P_; }

    void set_control_points(std::array<point_type, 4> const& P) { P_ = P; }

  private:
    std::array<point_type, 4> P_;
};

/**
 * @brief Curve formed of multiple bezier curves. The bezier curves should be connected.
 */
template <unsigned int Dimension, unsigned int Order>
class piecewise_bezier_curve_t
{
  public:
    using curve_segment_type = bezier_curve_t<Dimension, Order>;
    using point_type         = typename curve_segment_type::point_type;

    piecewise_bezier_curve_t() = default;
    piecewise_bezier_curve_t(std::vector<curve_segment_type> const& curve_segments)
        : curve_segments_(curve_segments)
    {
    }

    std::optional<point_type> eval(scalar_type t) const
    {
        assert(t >= 0.);
        index_type const i        = static_cast<index_type>(t);
        bool const is_t_in_bounds = i < curve_segments_.size();
        if (is_t_in_bounds)
        {
            t = t - static_cast<scalar_type>(i);
            return curve_segments_[i].eval(t);
        }
        return {};
    }

    void add_segment(curve_segment_type const& curve_segment)
    {
        if (!curve_segments_.empty())
        {
            auto const& previous_segment = curve_segments_.back();
            point_type const p0          = previous_segment.eval(1.);
            point_type const p1          = curve_segment.eval(0.);
            assert(p0.isApprox(p1, sbs::eps()));
        }
        curve_segments_.push_back(curve_segment);
    }

  private:
    std::vector<curve_segment_type> curve_segments_;
};

} // namespace geometry
} // namespace sbs

#endif // SBS_GEOMETRY_CURVES_H
