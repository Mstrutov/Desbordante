#pragma once

#include <optional>
#include <vector>

#include "algorithms/md/decision_boundary.h"
#include "algorithms/md/hymd/column_match_info.h"
#include "algorithms/md/hymd/indexes/pli_cluster.h"
#include "algorithms/md/hymd/indexes/records_info.h"
#include "algorithms/md/hymd/indexes/similarity_index.h"
#include "algorithms/md/hymd/indexes/similarity_matrix.h"
#include "algorithms/md/hymd/lattice/validation_info.h"
#include "algorithms/md/hymd/preprocessing/similarity.h"
#include "algorithms/md/hymd/preprocessing/similarity_measure/similarity_measure.h"
#include "model/index.h"
#include "util/worker_thread_pool.h"

namespace algos::hymd {

class SimilarityData {
public:
    using ColMatchesInfo = std::vector<
            std::tuple<std::unique_ptr<preprocessing::similarity_measure::SimilarityMeasure>,
                       model::Index, model::Index>>;

private:
    indexes::RecordsInfo const* const records_info_;
    bool const single_table_;

    std::vector<ColumnMatchInfo> const column_matches_sim_info_;
    std::vector<std::vector<model::md::DecisionBoundary>> const column_matches_lhs_bounds_;

    indexes::DictionaryCompressor const& GetLeftCompressor() const noexcept {
        return records_info_->GetLeftCompressor();
    }

    indexes::DictionaryCompressor const& GetRightCompressor() const noexcept {
        return records_info_->GetRightCompressor();
    }

public:
    SimilarityData(indexes::RecordsInfo* records_info,
                   std::vector<ColumnMatchInfo> column_matches_sim_info,
                   std::vector<std::vector<model::md::DecisionBoundary>>
                           column_matches_lhs_bounds) noexcept
        : records_info_(records_info),
          single_table_(records_info_->OneTableGiven()),
          column_matches_sim_info_(std::move(column_matches_sim_info)),
          column_matches_lhs_bounds_(std::move(column_matches_lhs_bounds)) {}

    static SimilarityData CreateFrom(indexes::RecordsInfo* records_info,
                                     ColMatchesInfo column_matches_info,
                                     util::WorkerThreadPool& pool);

    [[nodiscard]] std::size_t GetColumnMatchNumber() const noexcept {
        return column_matches_sim_info_.size();
    }

    std::vector<std::vector<model::md::DecisionBoundary>> const& GetLhsBounds() const noexcept {
        return column_matches_lhs_bounds_;
    }

    [[nodiscard]] std::vector<ColumnMatchInfo> const& GetColumnMatchesInfo() const noexcept {
        return column_matches_sim_info_;
    }

    [[nodiscard]] std::pair<model::Index, model::Index> GetColMatchIndices(
            model::Index index) const {
        auto const& [_, left_column_index, right_column_index] = column_matches_sim_info_[index];
        return {left_column_index, right_column_index};
    }

    [[nodiscard]] std::size_t GetLeftSize() const noexcept {
        return GetLeftCompressor().GetNumberOfRecords();
    }

    [[nodiscard]] std::unordered_set<PairComparisonResult> CompareAllWith(
            RecordIdentifier left_record_id) const;
    [[nodiscard]] PairComparisonResult CompareRecords(CompressedRecord const& left_record,
                                                      CompressedRecord const& right_record) const;

    [[nodiscard]] std::optional<model::md::DecisionBoundary> GetPreviousDecisionBound(
            model::md::DecisionBoundary lhs_bound, model::Index column_match_index) const;
};

}  // namespace algos::hymd
