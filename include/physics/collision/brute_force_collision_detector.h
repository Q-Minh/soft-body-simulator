#ifndef SBS_PHYSICS_COLLISION_BRUTE_FORCE_COLLISION_DETECTOR_H
#define SBS_PHYSICS_COLLISION_BRUTE_FORCE_COLLISION_DETECTOR_H

#include "common/primitive.h"
#include "physics/collision/collision_detector.h"

#include <vector>

namespace sbs {
namespace physics {

template <class ElementType, class ToCommonTetrahedronFunc>
class brute_force_collision_detector_t : public collision_detector_i<ElementType>
{
  public:
    brute_force_collision_detector_t(ToCommonTetrahedronFunc const& to_tetrahedron);

    void add_element(ElementType const& element);
    void add_element(ElementType&& element);
    template <class ElementIt>
    void add_elements(ElementIt first, ElementIt last);

    virtual std::list<ElementType> intersect(common::tetrahedron_t const& t) override;
    virtual std::list<ElementType> intersect(common::triangle_t const& t) override;

  private:
    ToCommonTetrahedronFunc to_tetrahedron_;
    std::vector<ElementType> elements_;
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_BRUTE_FORCE_COLLISION_DETECTOR_H