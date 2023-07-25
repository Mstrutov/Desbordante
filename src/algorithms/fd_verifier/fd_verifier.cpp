#include "algorithms/fd_verifier/fd_verifier.h"

#include <chrono>
#include <memory>
#include <stdexcept>

#include "util/config/equal_nulls/option.h"
#include "util/config/indices/option.h"
#include "util/config/indices/validate_index.h"
#include "util/config/names_and_descriptions.h"
#include "util/config/option_using.h"
#include "util/config/tabular_data/input_table/option.h"

namespace algos::fd_verifier {

FDVerifier::FDVerifier() : Algorithm({}) {
    RegisterOptions();
}

void FDVerifier::RegisterOptions() {
    DESBORDANTE_OPTION_USING;

    auto get_schema_cols = [this]() { return relation_->GetSchema()->GetNumColumns(); };
    auto check_rhs = [this](util::config::IndexType rhs_index) {
        util::config::ValidateIndex(rhs_index, relation_->GetSchema()->GetNumColumns());
    };

    RegisterInitialLoadOption(util::config::TableOpt(&input_table_));
    RegisterInitialLoadOption(util::config::EqualNullsOpt(&is_null_equal_null_));
    RegisterInitialExecOption(util::config::LhsIndicesOpt(&lhs_indices_, get_schema_cols));
    RegisterInitialExecOption(Option{&rhs_index_, kRhsIndex, kDRhsIndex}.SetValueCheck(check_rhs));
}

void FDVerifier::LoadDataInternal() {
    relation_ = ColumnLayoutRelationData::CreateFrom(*input_table_, is_null_equal_null_);
    input_table_->Reset();
    if (relation_->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty dataset: FD verifying is meaningless.");
    }
    typed_relation_ =
            model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_, is_null_equal_null_);
}

unsigned long long FDVerifier::ExecuteInternal() {
    auto start_time = std::chrono::system_clock::now();

    stats_calculator_ =
            std::make_unique<StatsCalculator>(relation_, typed_relation_, lhs_indices_, rhs_index_);

    VerifyFD();
    SortHighlightsByProportionDescending();
    stats_calculator_->PrintStatistics();

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    return elapsed_milliseconds.count();
}

void FDVerifier::VerifyFD() const {
    std::shared_ptr<util::PLI const> pli =
            relation_->GetColumnData(lhs_indices_[0]).GetPliOwnership();

    for (size_t i = 1; i < lhs_indices_.size(); ++i) {
        pli = pli->Intersect(relation_->GetColumnData(lhs_indices_[i]).GetPositionListIndex());
    }

    std::unique_ptr<util::PLI const> intersection_pli =
            pli->Intersect(relation_->GetColumnData(rhs_index_).GetPositionListIndex());

    if (pli->GetNumCluster() == intersection_pli->GetNumCluster()) {
        return;
    }

    auto& lhs_clusters = pli->GetIndex();
    stats_calculator_->CalculateStatistics(std::move(lhs_clusters));
}

void FDVerifier::SortHighlightsByProportionAscending() const {
    assert(stats_calculator_);
    stats_calculator_->SortHighlights(StatsCalculator::CompareHighlightsByProportionAscending());
}
void FDVerifier::SortHighlightsByProportionDescending() const {
    assert(stats_calculator_);
    stats_calculator_->SortHighlights(StatsCalculator::CompareHighlightsByProportionDescending());
}
void FDVerifier::SortHighlightsByNumAscending() const {
    assert(stats_calculator_);
    stats_calculator_->SortHighlights(StatsCalculator::CompareHighlightsByNumAscending());
}
void FDVerifier::SortHighlightsByNumDescending() const {
    assert(stats_calculator_);
    stats_calculator_->SortHighlights(StatsCalculator::CompareHighlightsByNumDescending());
}
void FDVerifier::SortHighlightsBySizeAscending() const {
    assert(stats_calculator_);
    stats_calculator_->SortHighlights(StatsCalculator::CompareHighlightsBySizeAscending());
}
void FDVerifier::SortHighlightsBySizeDescending() const {
    assert(stats_calculator_);
    stats_calculator_->SortHighlights(StatsCalculator::CompareHighlightsBySizeDescending());
}
void FDVerifier::SortHighlightsByLhsAscending() const {
    assert(stats_calculator_);
    stats_calculator_->SortHighlights(stats_calculator_->CompareHighlightsByLhsAscending());
}
void FDVerifier::SortHighlightsByLhsDescending() const {
    assert(stats_calculator_);
    stats_calculator_->SortHighlights(stats_calculator_->CompareHighlightsByLhsDescending());
}

}  // namespace algos::fd_verifier
