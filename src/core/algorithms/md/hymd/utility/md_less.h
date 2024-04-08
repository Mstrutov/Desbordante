#pragma once

#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/md.h"

namespace algos::hymd::utility {

using MdPair = std::pair<std::vector<model::md::DecisionBoundary>,
                         std::pair<model::Index, model::md::DecisionBoundary>>;

inline bool MdLessPairs(MdPair const& pair_left, MdPair const& pair_right) {
    auto const& [lhs_left, rhs_left] = pair_left;
    auto const& [lhs_right, rhs_right] = pair_right;
    auto const cardinality_left = std::count_if(lhs_left.begin(), lhs_left.end(),
                                                [](auto bound) { return bound != kLowestBound; });
    auto const cardinality_right = std::count_if(lhs_right.begin(), lhs_right.end(),
                                                 [](auto bound) { return bound != kLowestBound; });
    if (cardinality_left < cardinality_right) {
        return true;
    } else if (cardinality_left > cardinality_right) {
        return false;
    }
    auto comp = lhs_left <=> lhs_right;
    if (comp < 0) {
        return true;
    } else if (comp > 0) {
        return false;
    }
    auto [left_index, left_bound] = rhs_left;
    auto [right_index, right_bound] = rhs_right;
    if (left_bound > right_bound) {
        return true;
    } else if (left_bound < right_bound) {
        return false;
    }
    assert(left_index != right_index);
    return left_index < right_index;
}

inline bool MdLess(model::MD const& left, model::MD const& right) {
    return MdLessPairs(std::pair{left.GetLhsDecisionBounds(), left.GetRhs()},
                       std::pair{right.GetLhsDecisionBounds(), right.GetRhs()});
}
}  // namespace algos::hymd::utility
