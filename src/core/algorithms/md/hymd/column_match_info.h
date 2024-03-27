#pragma once

#include "algorithms/md/hymd/indexes/column_similarity_info.h"
#include "model/index.h"

namespace algos::hymd {

struct ColumnMatchInfo {
    indexes::ColumnMatchSimilarityInfo similarity_info;
    model::Index left_column_index;
    model::Index right_column_index;
};

}  // namespace algos::hymd
