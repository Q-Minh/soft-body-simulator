#include "common/mesh.h"

#include <Eigen/Core>
#include <numeric>

namespace sbs {
namespace common {

void shared_vertex_triangle_mesh_t::rescale(position_t const& boxmin, position_t const& boxmax)
{
    auto const min_reduce_op = [](position_t const& current_min, position_t const& p) {
        position_t new_min{current_min};
        if (p.x < current_min.x)
            new_min.x = p.x;
        if (p.y < current_min.y)
            new_min.y = p.y;
        if (p.z < current_min.z)
            new_min.z = p.z;

        return new_min;
    };

    auto const max_reduce_op = [](position_t const& current_max, position_t const& p) {
        position_t new_max{current_max};
        if (p.x > current_max.x)
            new_max.x = p.x;
        if (p.y > current_max.y)
            new_max.y = p.y;
        if (p.z > current_max.z)
            new_max.z = p.z;

        return new_max;
    };

    position_t const min_position = std::reduce(
        this->positions.begin(),
        this->positions.end(),
        this->positions.front(),
        min_reduce_op);

    position_t const max_position = std::reduce(
        this->positions.begin(),
        this->positions.end(),
        this->positions.front(),
        max_reduce_op);

    double const dx = max_position.x - min_position.x;
    double const dy = max_position.y - min_position.y;
    double const dz = max_position.z - min_position.z;

    double constexpr eps           = 1e-8;
    bool const dx_division_by_zero = std::abs(dx) < eps;
    bool const dy_division_by_zero = std::abs(dy) < eps;
    bool const dz_division_by_zero = std::abs(dz) < eps;

    auto const map_to_new_box = [=](position_t const& pos) {
        position_t rescaled_p{};
        rescaled_p.x = dx_division_by_zero ? pos.x : boxmin.x + (boxmax.x - boxmin.x) * (pos.x - min_position.x) / dx;
        rescaled_p.y = dy_division_by_zero ? pos.y : boxmin.y + (boxmax.y - boxmin.y) * (pos.y - min_position.y) / dy;
        rescaled_p.z = dz_division_by_zero ? pos.z : boxmin.z + (boxmax.z - boxmin.z) * (pos.z - min_position.z) / dz;
        return rescaled_p;
    };

    std::transform(
        this->positions.begin(),
        this->positions.end(),
        this->positions.begin(),
        [map_to_new_box](position_t const& pos) { return map_to_new_box(pos); });
}

void shared_vertex_tetrahedral_mesh_t::rescale(position_t const& boxmin, position_t const& boxmax)
{
    Eigen::Vector3d const box_minimum{boxmin.x, boxmin.y, boxmin.z};
    Eigen::Vector3d const box_maximum{boxmax.x, boxmax.y, boxmax.z};

    auto const min_reduce_op = [](position_t const& current_min, position_t const& p) {
        position_t new_min{current_min};
        if (p.x < current_min.x)
            new_min.x = p.x;
        if (p.y < current_min.y)
            new_min.y = p.y;
        if (p.z < current_min.z)
            new_min.z = p.z;

        return new_min;
    };

    auto const max_reduce_op = [](position_t const& current_max, position_t const& p) {
        position_t new_max{current_max};
        if (p.x > current_max.x)
            new_max.x = p.x;
        if (p.y > current_max.y)
            new_max.y = p.y;
        if (p.z > current_max.z)
            new_max.z = p.z;

        return new_max;
    };

    position_t const min_position = std::reduce(
        this->positions.begin(),
        this->positions.end(),
        this->positions.front(),
        min_reduce_op);

    position_t const max_position = std::reduce(
        this->positions.begin(),
        this->positions.end(),
        this->positions.front(),
        max_reduce_op);

    auto const map_to_new_box = [=](position_t const& pos) {
        Eigen::Vector3d const pmin{min_position.x, min_position.y, min_position.z};
        Eigen::Vector3d const pmax{max_position.x, max_position.y, max_position.z};
        Eigen::Vector3d const p{pos.x, pos.y, pos.z};

        Eigen::Vector3d const rescaled_p =
            box_minimum +
            Eigen::Vector3d{
                (box_maximum - box_minimum).array() * ((p - pmin).array() / (pmax - pmin).array())};

        return position_t{rescaled_p.x(), rescaled_p.y(), rescaled_p.z()};
    };

    std::transform(
        this->positions.begin(),
        this->positions.end(),
        this->positions.begin(),
        [map_to_new_box](position_t const& pos) { return map_to_new_box(pos); });
}

} // namespace common
} // namespace sbs