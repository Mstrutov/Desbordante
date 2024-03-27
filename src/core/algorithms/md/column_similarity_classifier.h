#pragma once

#include <cstddef>

#include "algorithms/md/decision_boundary.h"
#include "model/index.h"

namespace model::md {

class ColumnSimilarityClassifier {
private:
    Index column_match_index_;
    DecisionBoundary decision_boundary_;

public:
    ColumnSimilarityClassifier(Index column_match_index,
                               DecisionBoundary decision_boundary) noexcept
        : column_match_index_(column_match_index), decision_boundary_(decision_boundary) {}

    [[nodiscard]] Index GetColumnMatchIndex() const noexcept {
        return column_match_index_;
    }

    [[nodiscard]] DecisionBoundary GetDecisionBoundary() const noexcept {
        return decision_boundary_;
    }
};

}  // namespace model::md
