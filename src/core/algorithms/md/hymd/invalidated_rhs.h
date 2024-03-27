#pragma once

#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "model/index.h"

namespace algos::hymd {
struct InvalidatedRhs {
    model::Index index;
    model::md::DecisionBoundary old_bound;
    model::md::DecisionBoundary new_bound;
};

using InvalidatedRhss = std::vector<InvalidatedRhs>;
}  // namespace algos::hymd
