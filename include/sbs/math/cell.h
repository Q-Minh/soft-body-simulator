#ifndef SBS_MATH_CELL_H
#define SBS_MATH_CELL_H

#include "sbs/aliases.h"

#include <array>

namespace sbs {
namespace math {

template <unsigned int NumNodes, class BasisFunctionType>
class cell_t
{
  public:
    using node_count_value    = NumNodes;
    using size_type           = std::size_t;
    using basis_function_type = BasisFunctionType;

    // accessors
    size_type constexpr node_count() const { return nodes_.size(); }
    index_type node(index_type r) const { return nodes_[r]; }
    basis_function_type const& phi(index_type r) const { return phis_[r]; }
    basis_function_type& phi(index_type r) { return phis_[r]; }

    // mutators
    void set_phi(index_type r, basis_function_type const& phi) { phis_[r] = phi; }
    void set_node(index_type r, index_type i) { nodes_[r] = i; }

  private:
    std::array<index_type, NumNodes> nodes_;
    std::array<basis_function_type, NumNodes> phis_;
};

} // namespace math
} // namespace sbs

#endif // SBS_MATH_CELL_H