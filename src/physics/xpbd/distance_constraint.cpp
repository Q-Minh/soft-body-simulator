#include "physics/xpbd/distance_constraint.h"

#include "common/node.h"

#include <Eigen/Core>

namespace sbs {
namespace physics {
namespace xpbd {

void distance_constraint_t::project(
    common::scene_t& scene,
    scalar_type& lagrange_multiplier,
    scalar_type const dt) const
{
    auto const v1       = v1_.first;
    auto const v1_model = v1_.second;
    auto const v2       = v2_.first;
    auto const v2_model = v2_.second;

    //auto model1 =
    //    std::dynamic_pointer_cast<physics::node_t, common::node_t>(scene.objects[v1_model]);
    //auto model2 =
    //    std::dynamic_pointer_cast<physics::node_t, common::node_t>(scene.objects[v2_model]);

    //auto const& p1 = model1->get_positions()[v1];
    //auto const& p2 = model2->get_positions()[v2];

    //Eigen::Vector3d const p1_eigen{p1.x, p1.y, p1.z};
    //Eigen::Vector3d const p2_eigen{p2.x, p2.y, p2.z};

    //auto const diff = p2_eigen - p1_eigen;
    //auto const length = diff.norm();


}

} // namespace xpbd
} // namespace physics
} // namespace sbs