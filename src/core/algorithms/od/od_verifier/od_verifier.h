#pragma once

#include "algorithms/algorithm.h"
#include "algorithms/od/fastod/model/canonical_od.h"
#include "config/indices/type.h"
#include "model/table/column_layout_relation_data.h"

namespace algos::od_verifier {

class ODVerifier : public Algorithm {
private:
    using IndicesType = config::IndicesType;
    using IndexType = config::IndexType;
    using DataFrame = fastod::DataFrame;
    using PartitionCache = fastod::PartitionCache;
    using AscCanonicalOD = fastod::AscCanonicalOD;
    using DescCanonicalOD = fastod::DescCanonicalOD;
    using SimpleCanonicalOD = fastod::SimpleCanonicalOD;
    using AttributeSet = fastod::AttributeSet;

    // input data
    config::InputTable input_table_;
    config::EqNullsType is_null_equal_null_;
    IndicesType lhs_indices_;
    IndicesType rhs_indices_;
    IndicesType context_indices_;
    bool ascending_;

    // auxiliary data
    std::shared_ptr<ColumnLayoutRelationData> relation_;
    std::shared_ptr<DataFrame> data_;
    PartitionCache partition_cache_;

    // rows that vioalates ods
    std::vector<int> row_violate_ods_by_swap_;
    std::vector<int> row_violate_ods_by_split_;

    // load input data
    void RegisterOptions();
    void MakeExecuteOptsAvailable() override;
    void LoadDataInternal() override;

    // runs the algorithm and measures its time
    unsigned long long ExecuteInternal() override;

    // checks whether OD is violated and finds the rows where it is violated
    template <bool Ascending>
    void VerifyOD();

    // reset statistic of violations
    void ResetState() override {
        row_violate_ods_by_swap_.clear();
        row_violate_ods_by_split_.clear();
    }

public:
    // checks whether the OD has broken
    bool ODHolds() const {
        return row_violate_ods_by_swap_.empty() && row_violate_ods_by_split_.empty();
    }

    // base constructor
    ODVerifier();

    // Returns the number of rows that violate the OD by split
    size_t GetNumRowsViolateBySplit() const {
        return row_violate_ods_by_split_.size();
    }

    // Returns the number of rows that violate the OD by swap
    size_t GetNumRowsViolateBySwap() const {
        return row_violate_ods_by_swap_.size();
    }
};

}  // namespace algos::od_verifier
