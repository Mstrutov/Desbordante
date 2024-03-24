#pragma once

#include <vector>

#include "algorithms/algorithm.h"
#include "algorithms/ucc/ucc_verifier/ucc_stats_calculator.h"
#include "config/equal_nulls/type.h"
#include "config/indices/type.h"
#include "config/tabular_data/input_table_type.h"
#include "model/table/column_layout_relation_data.h"

namespace algos {

/* Algorithm used to verify that a set of columns is UCC and retrieving useful information in
 * case it is not */
class UCCVerifier : public Algorithm {
private:
    /* input options */
    config::IndicesType column_indices_;
    config::EqNullsType is_null_equal_null_{};

    config::InputTable input_table_;
    std::shared_ptr<ColumnLayoutRelationData> relation_;
    std::unique_ptr<UCCStatsCalculator> stats_calculator_;
    /* results of work */
    size_t num_rows_violating_ucc_ = 0;
    std::vector<model::PLI::Cluster> clusters_violating_ucc_;

    void VerifyUCC();
    void CalculateStatistics(std::deque<model::PLI::Cluster> const& clusters);
    void RegisterOptions();
    void LoadDataInternal() override;
    void MakeExecuteOptsAvailable() override;
    unsigned long long ExecuteInternal() override;
    std::shared_ptr<model::PLI const> CalculatePLI();

    void ResetState() override {
        if (stats_calculator_) {
            stats_calculator_->ResetState();
        }
    }

public:
    /* Returns true if UCC holds and false otherwise */
    bool UCCHolds() const {
        assert(stats_calculator_);
        return stats_calculator_->UCCHolds();
    }

    /* Returns the number of clusters where the UCC is violated, that is, the number of sets of rows
     * where each set consists of rows equal to each other in the specified columns */
    size_t GetNumClustersViolatingUCC() const {
        assert(stats_calculator_);
        return stats_calculator_->GetNumClustersViolatingUCC();
    }

    /* Returns the total number of table rows that violate the UCC */
    size_t GetNumRowsViolatingUCC() const {
        assert(stats_calculator_);
        return stats_calculator_->GetNumRowsViolatingUCC();
    }

    /* Returns clusters where the UCC is violated, that is, sets of rows where each set consists of
     * rows equal to each other in the specified columns */
    std::vector<model::PLI::Cluster> const& GetClustersViolatingUCC() const {
        assert(stats_calculator_);
        return stats_calculator_->GetClustersViolatingUCC();
    }

    double GetError() {
        assert(stats_calculator_);
        return stats_calculator_->GetAUCCError();
    }

    UCCVerifier();
};

}  // namespace algos
