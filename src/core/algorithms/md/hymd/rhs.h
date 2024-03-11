#pragma once

#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "model/index.h"

namespace algos::hymd {
struct Rhs {
    model::Index index;
    model::md::DecisionBoundary decision_boundary;
};

using Rhss = std::vector<Rhs>;
}  // namespace algos::hymd
