#include <cassert>
#include <sbs/common/mesh.h>
#include <sbs/physics/collision/bvh_model.h>
#include <sbs/physics/collision/contact.h>
#include <sbs/physics/collision/sdf_model.h>
#include <vector>

namespace sbs {
namespace physics {
namespace collision {

bool point_bvh_to_sdf_intersector_t::intersectVolume(Eigen::AlignedBox3d const& aabb) const
{
    Eigen::Vector3d const& min = aabb.min();
    Eigen::Vector3d const& max = aabb.max();

    auto const min_depth = sdf_.evaluate(min).first;
    auto const max_depth = sdf_.evaluate(max).first;

    bool const is_min_point_penetrating = min_depth < 0.;
    bool const is_max_point_penetrating = max_depth < 0.;

    return is_min_point_penetrating || is_max_point_penetrating;
}

bool point_bvh_to_sdf_intersector_t::intersectObject(index_type const vi) const
{
    Eigen::Vector3d const& point = surface_.vertex(vi).position;
    auto const [sdf, grad]       = sdf_.evaluate(point);

    bool const is_sdf_positive = sdf > 0.;
    if (is_sdf_positive)
        return false;

    scalar_type const penetration_depth   = sdf;
    Eigen::Vector3d const& contact_normal = grad;
    Eigen::Vector3d const contact_point   = point + -sdf * contact_normal;

    surface_mesh_particle_contact_t contact(
        contact_t::type_t::surface_particle_to_sdf,
        sdf_model_id_,
        bvh_model_id_,
        contact_point,
        contact_normal,
        vi);

    handler_.handle(contact);

    return false;
}

model_type_t point_bvh_model_t::model_type() const
{
    return model_type_t::bvh;
}

primitive_type_t point_bvh_model_t::primitive_type() const
{
    return primitive_type_t::point;
}

point_bvh_model_t::point_bvh_model_t(
    common::shared_vertex_surface_mesh_i const* surface)
    : surface_(surface)
{
    std::vector<index_type> vertices{};
    vertices.reserve(surface_->vertex_count());
    std::vector<Eigen::AlignedBox3d> aabbs{};
    aabbs.reserve(surface_->vertex_count());

    for (std::size_t i = 0u; i < surface_->vertex_count(); ++i)
    {
        Eigen::Vector3d const& p = surface_->vertex(i).position;
        Eigen::AlignedBox3d const aabb(p);

        vertices.push_back(static_cast<index_type>(i));
        aabbs.push_back(aabb);
    }

    bvh_.init(vertices.begin(), vertices.end(), aabbs.begin(), aabbs.end());

    /**
     * Balanced binary tree should have 2n-1 nodes in total.
     * There are n leaves, 1 root, and (n - 2) internal nodes.
     */
    auto const num_nodes_in_tree  = bvh_.getChildren().size();
    auto const num_internal_nodes = bvh_.getBoxes().size();
    auto const n                  = surface_->vertex_count();
    assert(num_nodes_in_tree == (2 * n - 1u));
    assert(num_internal_nodes == (n - 1u));
}

void point_bvh_model_t::collide(collision_model_t& other, contact_handler_t& handler)
{
    model_type_t const other_model_type = other.model_type();

    if (other_model_type == model_type_t::sdf)
    {
        sdf_model_t const& sdf_model = reinterpret_cast<sdf_model_t const&>(other);
        point_bvh_to_sdf_intersector_t
            intersector(*surface_, sdf_model, handler, this->id(), other.id());

        Eigen::BVIntersect(bvh_, intersector);
    }
}

void point_bvh_model_t::update(simulation_t const& simulation)
{
    auto const& boxes                = bvh_.getBoxes();
    auto const num_boxes             = boxes.size();
    auto const children_start_index  = num_boxes;
    std::vector<int> const& children = bvh_.getChildren();
    std::vector<index_type, Eigen::aligned_allocator<index_type>> const& objects =
        bvh_.getObjects();
    std::vector<int> bottom_up_traversal_queue{};
    bottom_up_traversal_queue.reserve(children.size() / 2u);
    // Update parents of leaf nodes and add them to the traversal queue for later
    for (std::size_t i = children_start_index; i < children.size(); i += 2)
    {
        int const parent_idx = static_cast<int>(children_start_index / 2u);
        assert(parent_idx < boxes.size());

        int const child1_idx = children[i];

        Eigen::AlignedBox3d& aabb = bvh_.getVolume(parent_idx);
        Eigen::Vector3d const& p1 = surface_->vertex(objects[child1_idx]).position;

        aabb.min() = p1;
        aabb.max() = p1;

        bool const has_second_child = i < children.size() - 1;
        if (has_second_child)
        {
            int const child2_idx      = children[i + 1];
            Eigen::Vector3d const& p2 = surface_->vertex(objects[child2_idx]).position;
            Eigen::Matrix<scalar_type, 3, 2> M{};
            M.col(0)   = p1;
            M.col(1)   = p2;
            aabb.min() = M.rowwise().minCoeff();
            aabb.max() = M.rowwise().maxCoeff();
        }

        bottom_up_traversal_queue.push_back(parent_idx);
    }

    // Traverse every level
    std::vector<int> bottom_up_traversal_queue_swap_buffer{};
    bottom_up_traversal_queue_swap_buffer.reserve(bottom_up_traversal_queue.size());
    // Until we are at the root node, we continue updating parent aabbs
    while (bottom_up_traversal_queue.size() > 1u)
    {
        bool const is_power_of_two = (bottom_up_traversal_queue.size() % 2 == 0);
        assert(is_power_of_two);
        for (std::size_t i = 0u; i < bottom_up_traversal_queue.size(); i += 2u)
        {
            int const child_idx              = bottom_up_traversal_queue[i];
            int const parent_idx             = child_idx / 2;
            Eigen::AlignedBox3d& child_aabb  = bvh_.getVolume(child_idx);
            Eigen::AlignedBox3d& parent_aabb = bvh_.getVolume(parent_idx);
            Eigen::Matrix<scalar_type, 3, 4> M{};
            M.col(0) = child_aabb.min();
            M.col(1) = parent_aabb.min();
            M.col(2) = child_aabb.max();
            M.col(3) = parent_aabb.max();

            parent_aabb.min() = M.rowwise().minCoeff();
            parent_aabb.max() = M.rowwise().maxCoeff();

            bottom_up_traversal_queue_swap_buffer.push_back(parent_idx);
        }
        bottom_up_traversal_queue.swap(bottom_up_traversal_queue_swap_buffer);
        bottom_up_traversal_queue_swap_buffer.clear();
    }
}

} // namespace collision
} // namespace physics
} // namespace sbs