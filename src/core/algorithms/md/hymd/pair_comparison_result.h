#pragma once

#include <vector>

#include "algorithms/md/hymd/preprocessing/similarity.h"
#include "algorithms/md/hymd/utility/vector_double_hash.h"

namespace algos::hymd {
using PairComparisonResult = std::vector<preprocessing::Similarity>;
}  // namespace algos::hymd
