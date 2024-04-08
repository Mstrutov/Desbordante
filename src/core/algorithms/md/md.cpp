#include "algorithms/md/md.h"

#include "model/table/column.h"

namespace model {

MD::MD(RelationalSchema const* left_schema, RelationalSchema const* right_schema,
       std::vector<md::ColumnMatch> column_matches,
       std::vector<md::LhsColumnSimilarityClassifier> lhs,
       md::ColumnSimilarityClassifier rhs) noexcept
    : left_schema_(left_schema),
      right_schema_(right_schema),
      column_matches_(std::move(column_matches)),
      lhs_(std::move(lhs)),
      rhs_(rhs) {}

std::string MD::ToStringFull() const {
    std::stringstream ss;
    ss << "[";
    for (auto const& classifier : lhs_) {
        // if (classifier.GetDecisionBoundary() == 0.0) continue;
        auto const& column_match = column_matches_[classifier.GetColumnMatchIndex()];
        ss << " ";
        ss << column_match.similarity_function_name << "(" << left_schema_->GetName() << ":"
           << left_schema_->GetColumn(column_match.left_col_index)->GetName() << ", "
           << right_schema_->GetName() << ":"
           << right_schema_->GetColumn(column_match.right_col_index)->GetName()
           << ")>=" << classifier.GetDecisionBoundary() << " ";
        auto disp = classifier.GetMaxDisprovedBound();
        if (disp.has_value()) {
            ss << "(>" << disp.value() << ") ";
        }
        ss << "|";
    }
    ss.seekp(-1, std::stringstream::cur);
    ss << "] -> ";
    auto const& classifier = rhs_;
    auto const& column_match = column_matches_[classifier.GetColumnMatchIndex()];
    ss << column_match.similarity_function_name << "(" << left_schema_->GetName() << ":"
       << left_schema_->GetColumn(column_match.left_col_index)->GetName() << ", "
       << right_schema_->GetName() << ":"
       << right_schema_->GetColumn(column_match.right_col_index)->GetName()
       << ")>=" << classifier.GetDecisionBoundary();
    return ss.str();
}

std::string MD::ToStringShort() const {
    std::stringstream ss;
    ss << "[";
    for (auto const& classifier : lhs_) {
        model::md::DecisionBoundary const decision_boundary = classifier.GetDecisionBoundary();
        if (decision_boundary != 0.0) {
            ss << decision_boundary;
        }
        ss << ",";
    }
    ss.seekp(-1, std::stringstream::cur);
    ss << "]->";
    ss << rhs_.GetColumnMatchIndex() << "@" << rhs_.GetDecisionBoundary();
    return ss.str();
}

std::vector<md::DecisionBoundary> MD::GetLhsDecisionBounds() const {
    std::vector<md::DecisionBoundary> decision_bounds;
    decision_bounds.reserve(lhs_.size());
    std::transform(
            lhs_.begin(), lhs_.end(), std::back_inserter(decision_bounds),
            [](md::LhsColumnSimilarityClassifier const& lhs) { return lhs.GetDecisionBoundary(); });
    return decision_bounds;
}

std::pair<Index, md::DecisionBoundary> MD::GetRhs() const noexcept {
    return {rhs_.GetColumnMatchIndex(), rhs_.GetDecisionBoundary()};
}

}  // namespace model
