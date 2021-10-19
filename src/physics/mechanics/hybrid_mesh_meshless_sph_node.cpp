#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_body.h>
#include <sbs/physics/mechanics/hybrid_mesh_meshless_sph_node.h>

namespace sbs {
namespace physics {
namespace mechanics {

hybrid_mesh_meshless_sph_node_t::hybrid_mesh_meshless_sph_node_t(
    index_type const i,
    functions::poly6_kernel_t const& kernel,
    hybrid_mesh_meshless_sph_body_t& body)
    : meshless_sph_node_t(i, kernel), body_(body), ti_(std::numeric_limits<index_type>::max())
{
}

void hybrid_mesh_meshless_sph_node_t::initialize(
    Eigen::Vector3d const& xi,
    std::vector<Eigen::Vector3d const*> const& pj,
    std::vector<index_type> const& neighbours,
    index_type const ti)
{
    meshless_sph_node_t::initialize(xi, pj, neighbours);
    ti_ = ti;
}

bool hybrid_mesh_meshless_sph_node_t::is_mixed_particle() const
{
    if (ti_ == std::numeric_limits<index_type>::max())
        return false;

    bool const is_in_boundary_tetrahedron = body_.is_boundary_mesh_tetrahedron(ti_);
    return is_in_boundary_tetrahedron;
}

hybrid_mesh_meshless_sph_body_t const& hybrid_mesh_meshless_sph_node_t::body() const
{
    return body_;
}

index_type hybrid_mesh_meshless_sph_node_t::ti() const
{
    return ti_;
}

} // namespace mechanics
} // namespace physics
} // namespace sbs