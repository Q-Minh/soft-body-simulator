#ifndef SBS_PHYSICS_COLLISION_COLLISION_DETECTOR_H
#define SBS_PHYSICS_COLLISION_COLLISION_DETECTOR_H

#include <list>

namespace sbs {

// forward declares
namespace common {

struct tetrahedron_t;
struct triangle_t;

} // namespace common

namespace physics {

template <class ElementType>
class collision_detector_i
{
  public:
    virtual std::list<ElementType> intersect(common::tetrahedron_t const& t) = 0;
    virtual std::list<ElementType> intersect(common::triangle_t const& t)    = 0;
};

} // namespace physics
} // namespace sbs

#endif // SBS_PHYSICS_COLLISION_COLLISION_DETECTOR_H