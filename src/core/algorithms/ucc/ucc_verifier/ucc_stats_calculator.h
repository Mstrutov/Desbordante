#pragma once

#include <vector>

#include "config/equal_nulls/type.h"
#include "config/indices/type.h"
#include "config/tabular_data/input_table_type.h"
#include "model/table/column_layout_relation_data.h"

namespace algos {
class UCCStatsCalculator {
  private:
    std::shared_ptr<ColumnLayoutRelationData>
        relation_;
    size_t num_rows_;
    double aucc_error_ = 0.0;
    /* results of work */
    size_t num_rows_violating_ucc_ = 0;
    std::vector<model::PLI::Cluster> clusters_violating_ucc_;

  public:
    UCCStatsCalculator(std::shared_ptr<ColumnLayoutRelationData> relation)
        : num_rows_(relation->GetNumRows()) {}
    UCCStatsCalculator(size_t num_rows) : num_rows_(num_rows) {}
    void ResetState() {
        num_rows_violating_ucc_ = 0;
        clusters_violating_ucc_.clear();
    }

    void CalculateStatistics(std::deque<model::PLI::Cluster> const &clusters) {
        // size_t num_rows = relation_->GetNumRows();

        unsigned long long num_pairs_combinations =
            static_cast<unsigned long long>(num_rows_);
        if (num_rows_ > 1) {
            num_pairs_combinations *= (num_rows_ - 1);
        }

        for (auto const &cluster : clusters) {
            num_rows_violating_ucc_ += cluster.size();
            clusters_violating_ucc_.push_back(cluster);
            aucc_error_ += static_cast<double>(cluster.size()) *
                           (cluster.size() - 1) / num_pairs_combinations;
        }
    }

    bool UCCHolds() const { return clusters_violating_ucc_.empty(); }

    /* Returns the number of clusters where the UCC is violated, that is, the
     * number of sets of rows where each set consists of rows equal to each
     * other in the specified columns */
    size_t GetNumClustersViolatingUCC() const {
        return clusters_violating_ucc_.size();
    }

    /* Returns the total number of table rows that violate the UCC */
    size_t GetNumRowsViolatingUCC() const { return num_rows_violating_ucc_; }

    /* Returns clusters where the UCC is violated, that is, sets of rows where
     * each set consists of rows equal to each other in the specified columns */
    std::vector<model::PLI::Cluster> const &GetClustersViolatingUCC() const {
        return clusters_violating_ucc_;
    }

    /* Returns error for aucc to hold*/
    double GetAUCCError() { return aucc_error_; }
};
} // namespace algos
