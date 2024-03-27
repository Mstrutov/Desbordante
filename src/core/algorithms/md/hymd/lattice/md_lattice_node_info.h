#pragma once

#include "algorithms/md/hymd/decision_boundary_vector.h"

namespace algos::hymd::lattice {

struct MdLatticeNodeInfo {
    DecisionBoundaryVector lhs_bounds;
    DecisionBoundaryVector* rhs_bounds;
};

}  // namespace algos::hymd::lattice
