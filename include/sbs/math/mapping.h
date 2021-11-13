#ifndef SBS_MATH_MAPPING_H
#define SBS_MATH_MAPPING_H

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace sbs {
namespace math {

struct tetrahedron_affine_mapping_t
{
    tetrahedron_affine_mapping_t(
        autodiff::Vector3dual X1,
        autodiff::Vector3dual X2,
        autodiff::Vector3dual X3,
        autodiff::Vector3dual X4,
        autodiff::Vector3dual x1,
        autodiff::Vector3dual x2,
        autodiff::Vector3dual x3,
        autodiff::Vector3dual x4)
        : D(), x1(x1)
    {
        D.col(0) = (x2 - x1);
        D.col(1) = (x3 - x1);
        D.col(2) = (x4 - x1);

        Dinv = D.inverse();
        det  = D.determinant();
    }

    autodiff::Vector3dual to_domain(autodiff::Vector3dual const& X) const { return x1 + D * X; }
    autodiff::Vector3dual to_ref(autodiff::Vector3dual const& x) const
    {
        autodiff::Vector3dual b = x - x1;
        return Dinv * b;
    }
    autodiff::Matrix3dual jacobian(autodiff::Vector3dual const& x) const { return D; }
    autodiff::dual determinant(autodiff::Vector3dual const& x) const { return det; }

    autodiff::Matrix3dual D;
    autodiff::Matrix3dual Dinv;
    autodiff::dual det;
    autodiff::Vector3dual x1;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_MAPPING_H