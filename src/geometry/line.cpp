#include "sbs/geometry/line.h"

#include <cmath>

namespace sbs {
namespace geometry {

std::tuple<std::vector<Eigen::Vector3d>, std::vector<index_type>, line_segment_t<3>> swept_surface(
    line_segment_t<3> const& line,
    std::function<std::optional<Eigen::Vector3d>(scalar_type)> const& parametric_curve,
    scalar_type t0,
    scalar_type dt,
    index_type n)
{
    std::optional<Eigen::Vector3d> const p0 = parametric_curve(t0);
    if (!p0.has_value())
        return std::make_tuple(std::vector<Eigen::Vector3d>{}, std::vector<index_type>{}, line);

    assert(line.p1().isApprox(*p0));

    std::pair<std::vector<Eigen::Vector3d>, std::vector<index_type>> surface{};

    Eigen::Vector3d const direction = line.direction();
    scalar_type const length        = line.length();
    line_segment_t<3> previous_line = line;
    for (auto i = 0u; i < n; ++i)
    {
        scalar_type const t                     = t0 + static_cast<scalar_type>(i + 1u) * dt;
        std::optional<Eigen::Vector3d> const Pi = parametric_curve(t);
        if (!Pi.has_value())
            break;

        line_segment_t<3> const next_line(*Pi, direction, length);

        auto [points, indices]  = previous_line.swept_surface(next_line);
        auto const index_offset = surface.first.size();

        std::copy(points.begin(), points.end(), std::back_inserter(surface.first));
        for (index_type const idx : indices)
        {
            surface.second.push_back(static_cast<index_type>(index_offset + idx));
        }

        previous_line = next_line;
    }

    return std::make_tuple(surface.first, surface.second, previous_line);
}

swept_line_segment_t::swept_line_segment_t(line_segment_t<3> const& starting_line_segment)
    : starting_line_segment_(starting_line_segment), t_(0.), dt_(0.1), vertices_(), triangles_()
{
}

void swept_line_segment_t::sweep(scalar_type delta_t)
{
    auto const n         = static_cast<index_type>(std::round(delta_t / dt_));
    scalar_type const t0 = t_;
    t_ += static_cast<scalar_type>(n) * dt_;

    auto const [points, indices, next_line] =
        swept_surface(starting_line_segment_, parametric_curve_, t0, dt_, n);

    starting_line_segment_ = next_line;
    std::vector<vertex_type> new_vertices{};
    new_vertices.reserve(points.size());
    for (auto const& p : points)
    {
        vertex_type v{};
        v.position = p;
        v.color    = color_;
        v.normal.setZero();

        new_vertices.push_back(v);
    }

    auto const index_offset = vertices_.size();
    for (auto i = 0u; i < indices.size(); i += 3u)
    {
        auto const v1 = indices[i];
        auto const v2 = indices[i + 1];
        auto const v3 = indices[i + 2];

        Eigen::Vector3d const p1 = points[v1];
        Eigen::Vector3d const p2 = points[v2];
        Eigen::Vector3d const p3 = points[v3];

        Eigen::Vector3d const normal = (p2 - p1).cross(p3 - p1);
        new_vertices[v1].normal += normal;
        new_vertices[v2].normal += normal;
        new_vertices[v3].normal += normal;

        triangle_type f(
            static_cast<index_type>(index_offset + v1),
            static_cast<index_type>(index_offset + v2),
            static_cast<index_type>(index_offset + v3));
        triangles_.push_back(f);
    }

    for (auto& v : new_vertices)
    {
        v.normal.normalize();
        vertices_.push_back(v);
    }
}

void swept_line_segment_t::prepare_vertices_for_rendering()
{
    std::vector<float> cpu_buffer{};
    cpu_buffer.reserve(vertices_.size() * 9u);
    for (auto const& v : vertices_)
    {
        cpu_buffer.push_back(static_cast<float>(v.position.x()));
        cpu_buffer.push_back(static_cast<float>(v.position.y()));
        cpu_buffer.push_back(static_cast<float>(v.position.z()));
        cpu_buffer.push_back(static_cast<float>(v.normal.x()));
        cpu_buffer.push_back(static_cast<float>(v.normal.y()));
        cpu_buffer.push_back(static_cast<float>(v.normal.z()));
        cpu_buffer.push_back(v.color.x());
        cpu_buffer.push_back(v.color.y());
        cpu_buffer.push_back(v.color.z());
    }
    this->transfer_vertices_for_rendering(std::move(cpu_buffer));
}

void swept_line_segment_t::prepare_indices_for_rendering()
{
    std::vector<unsigned int> index_buffer{};
    index_buffer.reserve(triangles_.size() * 3u);
    for (auto const& f : triangles_)
    {
        index_buffer.push_back(f.vertices[0]);
        index_buffer.push_back(f.vertices[1]);
        index_buffer.push_back(f.vertices[2]);
    }
    this->transfer_indices_for_rendering(std::move(index_buffer));
}

} // namespace geometry
} // namespace sbs