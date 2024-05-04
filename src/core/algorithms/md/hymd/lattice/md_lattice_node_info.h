#pragma once

#include <algorithm>

#include "algorithms/md/hymd/decision_boundary_vector.h"
#include "algorithms/md/hymd/lowest_bound.h"
#include "algorithms/md/hymd/md_lhs.h"

namespace algos::hymd::lattice {

struct MdLatticeNodeInfo {
    MdLhs lhs;
    DecisionBoundaryVector* rhs_bounds;

    void ZeroRhs() {
        std::fill(rhs_bounds->begin(), rhs_bounds->end(), kLowestBound);
    }
};

}  // namespace algos::hymd::lattice
