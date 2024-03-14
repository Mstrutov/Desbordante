#pragma once

#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "model/index.h"

namespace algos::hymd {
struct MdElement {
    model::Index index;
    model::md::DecisionBoundary decision_boundary;
};
}  // namespace algos::hymd
