#pragma once

#include <cstddef>

#include "algorithms/md/hymd/decision_boundary_vector.h"
#include "model/index.h"

namespace algos::hymd::utility {
inline model::Index GetFirstNonZeroIndex(DecisionBoundaryVector const& dec_bound_vector,
                                         model::Index index) {
    std::size_t const dec_bound_vec_size = dec_bound_vector.size();
    while (index < dec_bound_vec_size) {
        if (dec_bound_vector[index] != 0.0) return index;
        ++index;
    }
    return index;
}
}  // namespace algos::hymd::utility
