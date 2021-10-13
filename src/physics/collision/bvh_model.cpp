#include <cassert>
#include <sbs/common/mesh.h>
#include <sbs/physics/collision/bvh_model.h>
#include <sbs/physics/collision/contact.h>
#include <sbs/physics/collision/sdf_model.h>
#include <vector>

namespace sbs {
namespace physics {
namespace collision {

point_bvh_model_t::point_bvh_model_t() : kd_tree_type(0), surface_() {}

model_type_t point_bvh_model_t::model_type() const
{
    return model_type_t::bvh;
}

primitive_type_t point_bvh_model_t::primitive_type() const
{
    return primitive_type_t::point;
}

point_bvh_model_t::point_bvh_model_t(common::shared_vertex_surface_mesh_i const* surface)
    : Discregrid::KDTree<Discregrid::BoundingSphere>(surface->vertex_count()),
      surface_(surface)
{
    kd_tree_type::construct();
}

void point_bvh_model_t::collide(collision_model_t& other, contact_handler_t& handler)
{
    model_type_t const other_model_type = other.model_type();

    if (other_model_type == model_type_t::sdf)
    {
        sdf_model_t const& sdf_model = reinterpret_cast<sdf_model_t const&>(other);

        auto const is_sphere_colliding_with_sdf =
            [this, &sdf_model](unsigned int node_idx, unsigned int depth) -> bool {
            Discregrid::BoundingSphere const& s = this->hull(node_idx);
            auto const [signed_distance, grad]  = sdf_model.evaluate(s.x());
            bool const is_sphere_penetrating    = s.r() > std::abs(signed_distance);
            return is_sphere_penetrating;
        };

        auto const contact_callback =
            [this, &sdf_model, &handler](unsigned int node_idx, unsigned int depth) {
                kd_tree_type::Node const& node = this->node(node_idx);
                if (!node.isLeaf())
                    return;

                for (auto i = node.begin; i < node.begin + node.n; ++i)
                {
                    index_type const vi                = m_lst[i];
                    Eigen::Vector3d const& pi          = surface_->vertex(vi).position;
                    auto const [signed_distance, grad] = sdf_model.evaluate(pi);
                    bool const is_vertex_penetrating   = signed_distance < 0.;

                    if (!is_vertex_penetrating)
                        continue;

                    Eigen::Vector3d const contact_normal = grad.normalized();
                    Eigen::Vector3d const contact_point =
                        pi + std::abs(signed_distance) * contact_normal;

                    surface_mesh_particle_to_sdf_contact_t contact(
                        contact_t::type_t::surface_particle_to_sdf,
                        sdf_model.id(),
                        this->id(),
                        contact_point,
                        contact_normal,
                        vi);

                    handler.handle(contact);
                }
            };

        traverseBreadthFirst(is_sphere_colliding_with_sdf, contact_callback);
    }
}

void point_bvh_model_t::update(simulation_t const& simulation)
{
    kd_tree_type::update();
}

Eigen::Vector3d point_bvh_model_t::entityPosition(unsigned int i) const
{
    Eigen::Vector3d const position = surface_->vertex(i).position;
    return position;
}

void point_bvh_model_t::computeHull(
    unsigned int b,
    unsigned int n,
    Discregrid::BoundingSphere& hull) const
{
    auto vertices_of_sphere = std::vector<Eigen::Vector3d>(n);
    for (unsigned int i = b; i < n + b; ++i)
        vertices_of_sphere[i - b] = surface_->vertex(m_lst[i]).position;

    Discregrid::BoundingSphere const s(vertices_of_sphere);

    hull.x() = s.x();
    hull.r() = s.r();
}

} // namespace collision
} // namespace physics
} // namespace sbs