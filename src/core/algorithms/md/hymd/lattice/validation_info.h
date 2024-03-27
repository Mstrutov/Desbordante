#pragma once

#include <unordered_set>

#include <boost/dynamic_bitset.hpp>

#include "algorithms/md/hymd/lattice/md_lattice_node_info.h"
#include "model/index.h"

namespace algos::hymd::lattice {

struct ValidationInfo {
    MdLatticeNodeInfo* node_info;
    boost::dynamic_bitset<> rhs_indices;
};

}  // namespace algos::hymd::lattice
