#include <Eigen/Geometry>
#include <array>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <cassert>
#include <iostream>
#include <sbs/math/basis_functions.h>
#include <sbs/math/elasticity.h>
#include <sbs/math/interpolation.h>
#include <sbs/math/mapping.h>
#include <sbs/math/quadrature.h>

std::tuple<
    sbs::math::linear_hat_basis_function_op_t,
    sbs::math::linear_hat_basis_function_op_t,
    sbs::math::linear_hat_basis_function_op_t,
    sbs::math::linear_hat_basis_function_op_t>
get_basis_functions(
    autodiff::Vector3dual X1,
    autodiff::Vector3dual X2,
    autodiff::Vector3dual X3,
    autodiff::Vector3dual X4,
    autodiff::Vector3dual x1,
    autodiff::Vector3dual x2,
    autodiff::Vector3dual x3,
    autodiff::Vector3dual x4)
{
    autodiff::Matrix4dual A;
    unsigned int constexpr order = 1;
    A.col(0)                     = sbs::math::polynomial3d<order>(X1);
    A.col(1)                     = sbs::math::polynomial3d<order>(X2);
    A.col(2)                     = sbs::math::polynomial3d<order>(X3);
    A.col(3)                     = sbs::math::polynomial3d<order>(X4);

    autodiff::Matrix4dual Ainv;
    Ainv = A.inverse();

    sbs::math::linear_hat_basis_function_op_t phi1(Ainv.row(0u)), phi2(Ainv.row(1u)),
        phi3(Ainv.row(2u)), phi4(Ainv.row(3u));

    return {phi1, phi2, phi3, phi4};
}

int main()
{
    autodiff::Vector3dual X1, X2, X3, X4;
    X1 << 0., 0., 0.;
    X2 << 1., 0., 0.;
    X3 << 0., 1., 0.;
    X4 << 0., 0., 1.;

    autodiff::Matrix3dual S;
    S = 2. * autodiff::Matrix3dual::Identity();

    autodiff::Vector3dual x1 = S * X1;
    autodiff::Vector3dual x2 = S * X2;
    autodiff::Vector3dual x3 = S * X3;
    autodiff::Vector3dual x4 = S * X4;

    auto const [phi1, phi2, phi3, phi4] = get_basis_functions(X1, X2, X3, X4, x1, x2, x3, x4);

    double young_modulus = 1e6;
    double poisson_ratio = 0.35;

    using basis_function_type   = sbs::math::linear_hat_basis_function_op_t;
    using interpolation_op_type = sbs::math::interpolation_op_t<basis_function_type>;
    using deformation_gradient_op_type =
        sbs::math::deformation_gradient_op_t<interpolation_op_type>;
    using strain_op_type = sbs::math::strain_op_t<deformation_gradient_op_type>;

    interpolation_op_type interpolate_op({x1, x2, x3, x4}, {phi1, phi2, phi3, phi4});
    sbs::math::deformation_gradient_op_t<interpolation_op_type> deformation_gradient_op(
        interpolate_op);
    sbs::math::strain_op_t<deformation_gradient_op_type> strain_op(deformation_gradient_op);
    sbs::math::strain_energy_density_op_t<strain_op_type> strain_energy_density_op(
        strain_op,
        young_modulus,
        poisson_ratio);

    autodiff::Vector3dual refX1(0., 0., 0.);
    autodiff::Vector3dual refX2(1., 0., 0.);
    autodiff::Vector3dual refX3(0., 1., 0.);
    autodiff::Vector3dual refX4(0., 0., 1.);

    sbs::math::tetrahedron_affine_mapping_t mapping(X1, X2, X3, X4);
    sbs::math::tetrahedron_1point_constant_quadrature_rule_t<
        sbs::math::tetrahedron_affine_mapping_t>
        quadrature_rule(mapping);

    using autodiff::at;
    using autodiff::gradient;
    using autodiff::wrt;

    for (auto i = 0u; i < quadrature_rule.points.size(); ++i)
    {
        autodiff::Vector3dual x;
        autodiff::Matrix3dual F, E;

        auto const f = [&x, &F, &E, &strain_energy_density_op](
                           autodiff::Vector3dual Xi,
                           autodiff::dual wi) -> autodiff::dual {
            return wi * strain_energy_density_op(Xi, x, F, E);
        };

        autodiff::Vector3dual Xi = quadrature_rule.points[i];
        autodiff::dual wi        = quadrature_rule.weights[i];

        autodiff::dual wiPsi = f(Xi, wi);

        autodiff::Vector3dual f1 = autodiff::gradient(f, wrt(interpolate_op.uis[0]), at(Xi, wi));
        autodiff::Vector3dual f2 = autodiff::gradient(f, wrt(interpolate_op.uis[1]), at(Xi, wi));
        autodiff::Vector3dual f3 = autodiff::gradient(f, wrt(interpolate_op.uis[2]), at(Xi, wi));
        autodiff::Vector3dual f4 = autodiff::gradient(f, wrt(interpolate_op.uis[3]), at(Xi, wi));

        std::cout << "Xi\n" << Xi << "\n";
        std::cout << "wi=" << wi << "\n";
        std::cout << "Psi:\n" << (wiPsi.val / wi.val) << "\n";
        std::cout << "f1:\n" << f1 << "\n";
        std::cout << "f2:\n" << f2 << "\n";
        std::cout << "f3:\n" << f3 << "\n";
        std::cout << "f4:\n" << f4 << "\n";

        F = deformation_gradient_op(Xi, x);
        autodiff::Matrix3dual dPsidF;
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                autodiff::dual dPsidFij =
                    autodiff::derivative(strain_energy_density_op, wrt(F(i, j)), at(F, E));
                dPsidF(i, j) = dPsidFij;
            }
        }

        std::cout << "Piola stress:\n" << dPsidF << "\n";
    }

    return 0;
}
