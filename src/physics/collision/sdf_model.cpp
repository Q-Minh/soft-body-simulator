#include "..\..\..\include\sbs\physics\collision\sdf_model.h"

#include <sbs/physics/collision/bvh_model.h>
#include <sbs/physics/collision/contact.h>
#include <sbs/physics/collision/sdf_model.h>

namespace sbs {
namespace physics {
namespace collision {

sdf_model_t::sdf_model_t(
    Eigen::AlignedBox3d const& domain,
    std::array<unsigned int, 3u> const& resolution)
    : sdf_(domain, resolution)
{
}

sdf_model_t::sdf_model_t(Discregrid::CubicLagrangeDiscreteGrid const& sdf) : sdf_(sdf) {}

sdf_model_t::sdf_model_t(analytic_sdf_type const& analytic_sdf, Eigen::AlignedBox3d const& volume)
    : sdf_(Eigen::AlignedBox3d(), {2, 2, 2} /*dummy values*/), analytic_sdf_(analytic_sdf)
{
    this->volume() = volume;
}

model_type_t sdf_model_t::model_type() const
{
    return model_type_t::sdf;
}

primitive_type_t sdf_model_t::primitive_type() const
{
    return primitive_type_t::sdf;
}

void sdf_model_t::collide(collision_model_t& other, contact_handler_t& handler)
{
    model_type_t const other_model_type = other.model_type();
    if (other_model_type == model_type_t::bvh)
    {
        point_bvh_model_t& bvh = reinterpret_cast<point_bvh_model_t&>(other);
        bvh.collide(*this, handler);
    }
}

void sdf_model_t::update(simulation_t const& simulation)
{
    // Currently a no-op. However, for rigid bodies, we need to have a world-to-local transformation
    // to update.
}

sdf_model_t sdf_model_t::from_plane(
    Eigen::Hyperplane<scalar_type, 3> const& plane,
    Eigen::AlignedBox3d const& volume)
{
    auto const analytic_sdf =
        [plane](Eigen::Vector3d const& pi) -> std::pair<scalar_type, Eigen::Vector3d> {
        auto const sd   = plane.signedDistance(pi);
        auto const grad = plane.normal();
        return std::make_pair(sd, grad);
    };
    sdf_model_t sdf{analytic_sdf, volume};
    return sdf;
}

std::pair<scalar_type, Eigen::Vector3d> sdf_model_t::evaluate(Eigen::Vector3d const& p) const
{
    if (analytic_sdf_)
        return analytic_sdf_(p);

    unsigned int constexpr sdf_idx = 0u;
    Eigen::Vector3d grad{};
    double const signed_distance = sdf_.interpolate(sdf_idx, p, &grad);
    return {static_cast<scalar_type>(signed_distance), grad};
}

} // namespace collision
} // namespace physics
} // namespace sbs