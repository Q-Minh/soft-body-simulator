#ifndef SBS_MATH_QUADRATURE_H
#define SBS_MATH_QUADRATURE_H

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <vector>

namespace sbs {
namespace math {

template <class MappingType>
struct tetrahedron_1point_constant_quadrature_rule_t
{
    using mapping_type = MappingType;

    tetrahedron_1point_constant_quadrature_rule_t(mapping_type const& mapping) : points(), weights()
    {
        auto constexpr ref_volume = 1. / 6.;
        // reference tetrahedron barycenter 1/4 * ((0,0,0) + (1,0,0) + (0,1,0) + (0,0,1))
        autodiff::Vector3dual X(0.25, 0.25, 0.25);

        // get the quadrature point in domain space
        autodiff::Vector3dual x = mapping.to_domain(X);
        points.push_back(x);

        // get the quadrature weight as volume of the tetrahedron in domain space
        // multiplied by the volume change ratio induced by the mapping x: X->x
        autodiff::dual mapping_jacobian_det = mapping.determinant(x);
        autodiff::dual w                    = ref_volume * mapping_jacobian_det;
        weights.push_back(w);
    }

    std::vector<autodiff::Vector3dual> points;
    std::vector<autodiff::dual> weights;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_QUADRATURE_H