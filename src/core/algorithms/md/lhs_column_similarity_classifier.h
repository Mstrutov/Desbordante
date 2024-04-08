#pragma once

#include <optional>

#include "algorithms/md/column_similarity_classifier.h"
#include "algorithms/md/decision_boundary.h"
#include "model/index.h"

namespace model::md {

// HyMD search space consists of all matching dependencies with natural decision
// boundaries. However, they are generally not sufficient to describe all MDs
// that are valid on a pair of datasets. One may lower the decision boundary in
// any given column similarity classifier to a value higher than the previous
// natural decision boundary and still get a valid MD. This class stores that
// decision boundary. If the algorithm is configured to reduce the number of
// decision boundaries, then this value may not have a practical use but its
// name stays correct.
class LhsColumnSimilarityClassifier : public ColumnSimilarityClassifier {
private:
    // Missing value means the decision boundary does not affect the validity
    // of the whole MD, values of this column match can be as dissimilar as
    // possible with the whole MD staying valid.
    std::optional<DecisionBoundary> max_disproved_bound_;

public:
    LhsColumnSimilarityClassifier(std::optional<DecisionBoundary> max_disproved_bound,
                                  Index column_match_index,
                                  DecisionBoundary decision_boundary) noexcept
        : ColumnSimilarityClassifier(column_match_index, decision_boundary),
          max_disproved_bound_(std::move(max_disproved_bound)) {}

    [[nodiscard]] std::optional<DecisionBoundary> GetMaxDisprovedBound() const noexcept {
        return max_disproved_bound_;
    }
};

}  // namespace model::md
