#include "od_verifier.h"

#include "ascending_od/option.h"
#include "config/equal_nulls/option.h"
#include "config/indices/od_context.h"
#include "config/indices/option.h"
#include "config/tabular_data/input_table/option.h"
#include "partition.h"

namespace algos::od_verifier {

ODVerifier::ODVerifier() : Algorithm({}) {
    RegisterOptions();
    MakeOptionsAvailable({config::kTableOpt.GetName(), config::kEqualNullsOpt.GetName()});
}

void ODVerifier::RegisterOptions() {
    auto get_schema_cols = [this]() { return relation_->GetSchema()->GetNumColumns(); };

    RegisterOption(config::kTableOpt(&input_table_));
    RegisterOption(config::kEqualNullsOpt(&is_null_equal_null_));
    RegisterOption(config::kLhsIndicesOpt(&lhs_indices_, get_schema_cols));
    RegisterOption(config::kRhsIndicesOpt(&rhs_indices_, get_schema_cols));
    RegisterOption(config::kODContextOpt(&context_indices_));
    RegisterOption(config::kAscendingODOpt(&ascending_));
}

void ODVerifier::MakeExecuteOptsAvailable() {
    MakeOptionsAvailable({config::kLhsIndicesOpt.GetName(), config::kRhsIndicesOpt.GetName(),
                          config::kODContextOpt.GetName(), config::kAscendingODOpt.GetName()});
}

void ODVerifier::LoadDataInternal() {
    relation_ = ColumnLayoutRelationData::CreateFrom(*input_table_, is_null_equal_null_);

    if (relation_->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty dataset: OD verifying is meaningless.");
    }
    input_table_->Reset();
    data_ = std::make_shared<DataFrame>(DataFrame::FromInputTable(input_table_));
    if (data_->GetColumnCount() == 0)
        throw std::runtime_error("Got an empty dataset: OD verifying is meaningless.");
}

unsigned long long ODVerifier::ExecuteInternal() {
    auto start_time = std::chrono::system_clock::now();
    if (ascending_)
        VerifyOD<true>();
    else
        VerifyOD<false>();
    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    return elapsed_milliseconds.count();
}

template <bool Ascending>
void ODVerifier::VerifyOD() {
    AttributeSet context;

    for (auto column : context_indices_) context.Set(column);

    fastod::ComplexStrippedPartition stripped_partition_swap(
            (partition_cache_.GetStrippedPartition(context, data_)));

    if (stripped_partition_swap.Swap<Ascending>(lhs_indices_[0], rhs_indices_[0])) {
        Partition part{stripped_partition_swap};
        std::vector<std::pair<int, int>> violates(
                part.FindViolationsBySwap<Ascending>(lhs_indices_[0], rhs_indices_[0]));

        for (auto position_violate : violates)
            row_violate_ods_by_swap_.push_back(position_violate.second + 1);
    }

    context.Set(lhs_indices_[0]);
    fastod::ComplexStrippedPartition stripped_partition_split(
            partition_cache_.GetStrippedPartition(context, data_));

    if (stripped_partition_split.Split(rhs_indices_[0])) {
        Partition part{stripped_partition_split};
        std::vector<std::pair<int, int>> violates(part.FindViolationsBySplit(rhs_indices_[0]));

        for (auto position_violate : violates)
            row_violate_ods_by_split_.push_back(position_violate.second + 1);
    }
    std::sort(row_violate_ods_by_split_.begin(), row_violate_ods_by_split_.end());
    std::sort(row_violate_ods_by_swap_.begin(), row_violate_ods_by_swap_.end());
}

}  // namespace algos::od_verifier
