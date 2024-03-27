#pragma once

#include <cstddef>
#include <functional>

#include "algorithms/md/decision_boundary.h"
#include "model/index.h"

namespace algos::hymd::lattice {
using SingleLevelFunc = std::function<std::size_t(model::md::DecisionBoundary, model::Index)>;
}  // namespace algos::hymd::lattice
