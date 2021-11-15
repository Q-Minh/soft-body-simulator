#include <Eigen/Geometry>
#include <array>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <cassert>
#include <iostream>
#include <sbs/geometry/get_simple_bar_model.h>
#include <sbs/geometry/grid.h>
#include <sbs/geometry/tetrahedral_domain.h>
#include <sbs/math/basis_functions.h>
#include <sbs/math/elasticity.h>
#include <sbs/math/interpolation.h>
#include <sbs/math/mapping.h>
#include <sbs/math/mls.h>
#include <sbs/math/quadrature.h>
#include <sbs/physics/mechanics/efg_tetrahedral_meshless_model.h>

int main()
{
    using kernel_function_type   = sbs::math::quartic_spline_kernel_t;
    unsigned int constexpr order = 1u;
    using meshless_model_type =
        sbs::physics::mechanics::efg_tetrahedral_meshless_model_t<kernel_function_type, order>;
    using basis_function_type         = typename meshless_model_type::basis_function_type;
    using interpolation_function_type = typename meshless_model_type::interpolation_function_type;

    sbs::scalar_type const hmultiplier             = 1.1;
    sbs::common::geometry_t const geometry         = sbs::geometry::get_simple_bar_model(2, 2, 2);
    std::vector<Eigen::Vector3d> const tet_points  = sbs::common::to_points(geometry);
    std::vector<sbs::index_type> const tet_indices = sbs::common::to_indices(geometry);
    sbs::geometry::tetrahedral_domain_t const domain(tet_points, tet_indices);
    Eigen::Vector3i const resolution{3, 3, 3};
    sbs::geometry::grid_t const grid(domain, resolution);

    meshless_model_type meshless_model(domain, grid, hmultiplier);
    auto const& topology         = domain.topology();
    auto const tetrahedron_count = topology.tetrahedron_count();
    for (sbs::index_type ti = 0u; ti < tetrahedron_count; ++ti)
    {
        sbs::topology::tetrahedron_t const& t = topology.tetrahedron(ti);
        Eigen::Vector3d const& p1             = domain.position(t.v1());
        Eigen::Vector3d const& p2             = domain.position(t.v2());
        Eigen::Vector3d const& p3             = domain.position(t.v3());
        Eigen::Vector3d const& p4             = domain.position(t.v4());

        // Insert integration point at the barycenter of the integration domain's tetrahedra

        sbs::math::tetrahedron_affine_mapping_t const mapping(p1, p2, p3, p4);
        sbs::math::tetrahedron_1point_constant_quadrature_rule_t<decltype(mapping)> quadrature_rule{
            mapping};

        for (auto i = 0u; i < quadrature_rule.points.size(); ++i)
        {
            auto const Xi = quadrature_rule.points[i];
            auto const wi = quadrature_rule.weights[i];

            Eigen::Vector3d const integration_point = Xi.cast<sbs::scalar_type>();
            sbs::index_type const integration_point_idx =
                static_cast<sbs::index_type>(meshless_model.integration_point_count());
            bool const was_integration_point_added =
                meshless_model.add_integration_point(integration_point);
        }
    }
    for (auto i = 0u; i < meshless_model.integration_point_count(); ++i)
    {
        Eigen::Vector3d const& integration_point = meshless_model.integration_point(i);
        interpolation_function_type const& interpolate =
            meshless_model.interpolation_field_from_integration_point(i);
        autodiff::Vector3dual const x = interpolate(integration_point);

        sbs::math::mls_deformation_gradient_op_t<interpolation_function_type>
            deformation_gradient_function(interpolate);
        autodiff::Matrix3dual F = deformation_gradient_function(integration_point);

        std::cout << "Integration point " << i << ":\n" << integration_point << "\n";
        std::cout << "x" << i << ":\n" << x << "\n";
        //std::cout << "F" << i << ":\n" << F << "\n\n";
    }

    return 0;
}
