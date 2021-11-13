#ifndef SBS_ALIASES_H
#define SBS_ALIASES_H

#include <cstdint>

namespace sbs {

using index_type  = std::uint32_t;
using scalar_type = double;

scalar_type constexpr eps()
{
    return 1e-15;
}

} // namespace sbs

#endif // SBS_ALIASES_H