#pragma once

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "algorithms/md/hymd/preprocessing/similarity.h"
#include "algorithms/md/hymd/table_identifiers.h"

namespace algos::hymd::indexes {
using RecSet = std::unordered_set<RecordIdentifier>;
using MatchingRecsMapping = std::map<preprocessing::Similarity, RecSet>;
using SimilarityIndex = std::vector<MatchingRecsMapping>;
}  // namespace algos::hymd::indexes
