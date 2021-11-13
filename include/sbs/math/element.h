#ifndef SBS_MATH_ELEMENT_H
#define SBS_MATH_ELEMENT_H

#include "mapping.h"

#include <array>
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

namespace sbs {
namespace math {

template <unsigned int i>
constexpr autodiff::Vector3dual reference_tetrahedron_position()
{
    static_assert(i >= 0 && i < 4);
    if constexpr (i == 0)
    {
        return autodiff::Vector3dual(0., 0., 0.);
    }
    if constexpr (i == 1)
    {
        return autodiff::Vector3dual(1., 0., 0.);
    }
    if constexpr (i == 2)
    {
        return autodiff::Vector3dual(0., 1., 0.);
    }
    if constexpr (i == 3)
    {
        return autodiff::Vector3dual(0., 0., 1.);
    }
}

std::array<autodiff::Vector3dual, 4> reference_tetrahedron()
{
    return {
        reference_tetrahedron_position<0>(),
        reference_tetrahedron_position<1>(),
        reference_tetrahedron_position<2>(),
        reference_tetrahedron_position<3>()};
}

class tetrahedral_element_t
{
  public:
    tetrahedral_element_t(
        autodiff::Vector3dual const& X1,
        autodiff::Vector3dual const& X2,
        autodiff::Vector3dual const& X3,
        autodiff::Vector3dual const& X4)
        : mapping_(
              reference_tetrahedron_position<0>(),
              reference_tetrahedron_position<1>(),
              reference_tetrahedron_position<2>(),
              reference_tetrahedron_position<3>(),
              X1,
              X2,
              X3,
              X4)
    {
    }

    tetrahedron_affine_mapping_t const& mapping() const { return mapping_; }

  private:
    tetrahedron_affine_mapping_t mapping_;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_ELEMENT_H