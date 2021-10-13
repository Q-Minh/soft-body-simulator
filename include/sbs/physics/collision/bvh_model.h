#ifndef SBS_PHYSICS_COLLISION_BVH_MODEL_H
#define SBS_PHYSICS_COLLISION_BVH_MODEL_H

#include <Discregrid/acceleration/bounding_sphere.hpp>
#include <Discregrid/acceleration/kd_tree.hpp>
#include <sbs/aliases.h>
#include <sbs/physics/collision/collision_model.h>

namespace sbs {
namespace common {

class shared_vertex_surface_mesh_i;
class contact_handler_t;

} // namespace common

namespace physics {

class simulation_t;

namespace collision {

class sdf_model_t;

class point_bvh_to_sdf_intersector_t
{
  public:
    point_bvh_to_sdf_intersector_t(
        common::shared_vertex_surface_mesh_i const& surface,
        sdf_model_t const& sdf,
        contact_handler_t& handler,
        index_type const bvh_model_id,
        index_type const sdf_model_id)
        : surface_(surface),
          sdf_(sdf),
          handler_(handler),
          bvh_model_id_(bvh_model_id),
          sdf_model_id_(sdf_model_id)
    {
    }

    bool intersectVolume(Eigen::AlignedBox3d const& aabb) const;
    bool intersectObject(index_type const vi) const;

  private:
    common::shared_vertex_surface_mesh_i const& surface_;
    sdf_model_t const& sdf_;
    contact_handler_t& handler_;
    index_type bvh_model_id_;
    index_type sdf_model_id_;
};

class point_bvh_model_t : public collision_model_t,
                          public Discregrid::KDTree<Discregrid::BoundingSphere>
{
  public:
    point_bvh_model_t();
    point_bvh_model_t(common::shared_vertex_surface_mesh_i const* surface);

    point_bvh_model_t(point_bvh_model_t const& other) = default;
    point_bvh_model_t(point_bvh_model_t&& other)      = default;
    point_bvh_model_t& operator=(point_bvh_model_t const& other) = default;
    point_bvh_model_t& operator=(point_bvh_model_t&& other) = default;

    virtual model_type_t model_type() const override;
    virtual primitive_type_t primitive_type() const override;
    virtual void collide(collision_model_t& other, contact_handler_t& handler) override;
    virtual void update(simulation_t const& simulation) override;

  protected:
    using kd_tree_type = Discregrid::KDTree<Discregrid::BoundingSphere>;

    virtual Eigen::Vector3d entityPosition(unsigned int i) const override final;
    virtual void computeHull(unsigned int b, unsigned int n, Discregrid::BoundingSphere& hull)
        const override final;

  private:
    common::shared_vertex_surface_mesh_i const* surface_;
    contact_handler_t* handler_;
};

} // namespace collision
} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_BVH_MODEL_H