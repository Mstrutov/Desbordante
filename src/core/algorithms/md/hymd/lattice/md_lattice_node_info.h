#pragma once

#include "algorithms/md/hymd/decision_boundary_vector.h"
#include "algorithms/md/hymd/md_lhs.h"

namespace algos::hymd::lattice {

struct MdLatticeNodeInfo {
    MdLhs lhs;
    DecisionBoundaryVector* rhs_bounds;
};

}  // namespace algos::hymd::lattice
