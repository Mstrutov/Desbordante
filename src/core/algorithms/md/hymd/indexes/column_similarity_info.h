#pragma once

#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/indexes/similarity_index.h"
#include "algorithms/md/hymd/indexes/similarity_matrix.h"
#include "algorithms/md/hymd/preprocessing/similarity.h"

namespace algos::hymd::indexes {
struct ColumnMatchSimilarityInfo {
    preprocessing::Similarity lowest_similarity;
    SimilarityMatrix similarity_matrix;
    SimilarityIndex similarity_index;
    // TODO: add slim similarity index
};

struct SimilarityMeasureOutput {
    std::vector<model::md::DecisionBoundary> lhs_bounds;
    ColumnMatchSimilarityInfo indexes;
};
}  // namespace algos::hymd::indexes
