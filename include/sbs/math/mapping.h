#ifndef SBS_MATH_MAPPING_H
#define SBS_MATH_MAPPING_H

#include "sbs/aliases.h"

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace sbs {
namespace math {

class identity_mapping_t
{
  public:
    autodiff::Vector3dual to_domain(autodiff::Vector3dual const& X) const { return X; }
    autodiff::Vector3dual to_ref(autodiff::Vector3dual const& x) const { return x; }
    autodiff::Matrix3dual jacobian(autodiff::Vector3dual const& x) const
    {
        return autodiff::Matrix3dual::Identity();
    }
    autodiff::dual determinant(autodiff::Vector3dual const& x) const { return 1.; }
};

class tetrahedron_affine_mapping_t
{
  public:
    tetrahedron_affine_mapping_t(
        autodiff::Vector3dual x1,
        autodiff::Vector3dual x2,
        autodiff::Vector3dual x3,
        autodiff::Vector3dual x4)
        : D_(), x1_(x1)
    {
        D_.col(0) = (x2 - x1);
        D_.col(1) = (x3 - x1);
        D_.col(2) = (x4 - x1);

        Dinv_ = D_.inverse();
        det_  = D_.determinant();
    }

    autodiff::Vector3dual to_domain(autodiff::Vector3dual const& X) const { return x1_ + D_ * X; }
    autodiff::Vector3dual to_ref(autodiff::Vector3dual const& x) const
    {
        autodiff::Vector3dual b = x - x1_;
        return Dinv_ * b;
    }
    autodiff::Matrix3dual jacobian(autodiff::Vector3dual const& x) const { return D_; }
    autodiff::dual determinant(autodiff::Vector3dual const& x) const { return det_; }

  private:
    autodiff::Matrix3dual D_;
    autodiff::Matrix3dual Dinv_;
    autodiff::dual det_;
    autodiff::Vector3dual x1_;
};

class tetrahedron_barycentric_mapping_t
{
  public:
    tetrahedron_barycentric_mapping_t(
        autodiff::Vector3dual const& X1,
        autodiff::Vector3dual const& X2,
        autodiff::Vector3dual const& X3,
        autodiff::Vector3dual const& X4)
        : A_(), Ainv_(), det_()
    {
        A_.row(0).setOnes();
        A_.block(1, 0, 3, 1) = X1.cast<scalar_type>();
        A_.block(1, 1, 3, 1) = X2.cast<scalar_type>();
        A_.block(1, 2, 3, 1) = X3.cast<scalar_type>();
        A_.block(1, 3, 3, 1) = X4.cast<scalar_type>();

        Ainv_ = A_.inverse();
        det_  = A_.determinant();
    }

    autodiff::Vector4dual barycentric_coordinates(autodiff::Vector3dual X) const
    {
        autodiff::Vector4dual P1{1., X.x(), X.y(), X.z()};
        autodiff::Vector4dual bc = Ainv_ * P1;
        return bc;
    }
    autodiff::Vector3dual coordinates(autodiff::Vector4dual bc) const
    {
        autodiff::Vector4dual P1 = A_ * bc;
        return autodiff::Vector3dual(P1(1), P1(2), P1(3));
    }
    autodiff::dual determinant() const { return det_; }

  private:
    Eigen::Matrix4d A_;
    Eigen::Matrix4d Ainv_;
    autodiff::dual det_;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_MAPPING_H