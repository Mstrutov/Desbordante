#pragma once

#include <vector>

#include <boost/unordered/unordered_flat_map.hpp>

#include "algorithms/md/hymd/preprocessing/similarity.h"
#include "algorithms/md/hymd/table_identifiers.h"

namespace algos::hymd::indexes {
using SimilarityMatrixRow =
        boost::unordered::unordered_flat_map<ValueIdentifier, preprocessing::Similarity>;
using SimilarityMatrix = std::vector<SimilarityMatrixRow>;
}  // namespace algos::hymd::indexes
