#include "cords.h"

#include <chrono>
#include <string>
#include <utility>
#include <vector>

#include <boost/math/distributions/chi_squared.hpp>

#include "config/equal_nulls/option.h"
#include "config/names_and_descriptions.h"
#include "config/option.h"
#include "config/tabular_data/input_table/option.h"
#include "frequency_handler.h"
#include "model/table/column_index.h"
#include "model/table/typed_column_data.h"
#include "sample.h"

namespace algos {
Cords::Cords() : FDAlgorithm({kFirstPhaseName, kSecondPhaseName}) {
    RegisterOptions();
    MakeOptionsAvailable({config::kTableOpt.GetName(), config::kEqualNullsOpt.GetName()});
}

void Cords::RegisterOptions() {
    using namespace config::names;
    using namespace config::descriptions;
    using config::Option;
    auto check_param = [](long double param) {
        if (param < 0 || param > 1) throw config::ConfigurationError("Parameter out of range");
    };
    auto check_max_false_positive = [](long double param) {
        if (param < 0 || param > 0.39)
            throw config::ConfigurationError(
                    "Maximum probability of a false-positive result is out of range");
    };
    auto check_positive = [](size_t param) {
        if (param == 0) throw config::ConfigurationError("Parameter out of range");
    };
    auto check_delta = [this](long double param) {
        if (param < 0 || param > 1) throw config::ConfigurationError("delta out of range");
        if (param < minimum_cardinality_) {
            throw config::ConfigurationError("delta must be less than minimum_cardinality_");
        }
    };
    RegisterOption(config::kTableOpt(&input_table_));
    RegisterOption(config::kEqualNullsOpt(&is_null_equal_null_));

    RegisterOption(Option{&only_sfd_, kOnlySFD, kDOnlySFD});
    RegisterOption(Option{&minimum_cardinality_, kMinCard, kDMinCard}.SetValueCheck(check_param));
    RegisterOption(
            Option{&max_diff_vals_proportion_, kMaxDiffValsProportion, kDMaxDiffValsProportion}
                    .SetValueCheck(check_param));
    RegisterOption(
            Option{&min_sfd_strength_measure_, kMinSFDStrengthMeasure, kDMinSFDStrengthMeasure}
                    .SetValueCheck(check_param));
    RegisterOption(
            Option{&min_skew_threshold_, kMinSkewThreshold, kDMinSkewThreshold}.SetValueCheck(
                    check_param));
    RegisterOption(Option{&min_structural_zeroes_proportion_, kMinStructuralZeroesAmount,
                          kDMinStructuralZeroesAmount}
                           .SetValueCheck(check_param));
    RegisterOption(Option{&max_false_positive_probability_, kMaxFalsePositiveProbability,
                          kDMaxFalsePositiveProbability}
                           .SetValueCheck(check_max_false_positive));

    RegisterOption(Option{&delta_, kDelta, kDDelta}.SetValueCheck(check_delta));
    RegisterOption(
            Option{&max_amount_of_categories_, kMaxAmountOfCategories, kDMaxAmountOfCategories}
                    .SetValueCheck(check_positive));
}

void Cords::MakeExecuteOptsAvailableFDInternal() {
    using namespace config::names;
    MakeOptionsAvailable({kOnlySFD, kMinCard, kMaxDiffValsProportion, kMinSFDStrengthMeasure,
                          kMinSkewThreshold, kMinStructuralZeroesAmount,
                          kMaxFalsePositiveProbability, kDelta, kMaxAmountOfCategories});
}

void Cords::ResetStateFd() {
    is_skewed_.clear();
    domains_.clear();
    soft_keys_.clear();
    trivial_columns_.clear();
    correlations_collection_.Clear();
}

void Cords::LoadDataInternal() {
    typed_relation_ =
            model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_, is_null_equal_null_);
}

bool Cords::DetectAndRegisterSFD(Sample const &smp) {
    if (smp.GetConcatCardinality() <= max_diff_vals_proportion_ * smp.GetRowIndices().size()) {
        if (smp.GetLhsCardinality() >=
            (1 - min_sfd_strength_measure_) * smp.GetConcatCardinality()) {
            RegisterFd(smp.GetLhsVertical(), smp.GetRhsColumn());
            return true;
        } else if (smp.GetRhsCardinality() >=
                   (1 - min_sfd_strength_measure_) * smp.GetConcatCardinality()) {
            RegisterFd(smp.GetRhsVertical(), smp.GetLhsColumn());
            return true;
        }
    }
    return false;
}

void Cords::SkewHandling(model::ColumnIndex col_i, model::ColumnIndex col_k,
                         std::vector<model::TypedColumnData> const &data, Sample &smp,
                         FrequencyHandler const &handler) {
    for (model::ColumnIndex col_ind : {col_i, col_k}) {
        if (handler.GetColumnFrequencySum(col_ind) >=
            (1 - min_skew_threshold_) * data[col_ind].GetNumRows()) {
            is_skewed_[col_ind] = true;
            domains_[col_ind] = handler.ColumnFrequencyMapSize(col_ind);
            smp.Filter(handler, data, col_ind);
        } else {
            domains_[col_ind] =
                    std::min(handler.GetColumnCardinality(col_ind), max_amount_of_categories_);
        }
    }
}

size_t Cords::Category(model::ColumnIndex col_ind, std::string const &val, size_t domain, bool skew,
                       FrequencyHandler const &handler) const {
    if (skew) {
        return handler.GetValueOrdinalNumberAtColumn(val, col_ind);
    }
    return std::hash<std::string>{}(val) % domain;
}

void Cords::FillContingencyTable(Sample &smp, std::vector<model::TypedColumnData> const &data,
                                 model::ColumnIndex col_i, model::ColumnIndex col_k,
                                 std::vector<std::vector<long double>> &n_i_j,
                                 std::vector<long double> &n_i, std::vector<long double> &n_j,
                                 FrequencyHandler const &handler) const {
    for (size_t row_ind : smp.GetRowIndices()) {
        size_t i = Category(col_i, data[col_i].GetDataAsString(row_ind), domains_[col_i],
                            is_skewed_[col_i], handler);
        size_t j = Category(col_k, data[col_k].GetDataAsString(row_ind), domains_[col_k],
                            is_skewed_[col_k], handler);
        n_i_j[i][j]++;
        n_i[i]++;
        n_j[j]++;
    }
}

bool Cords::TooMuchStructuralZeros(std::vector<std::vector<long double>> const &n_i_j,
                                   model::ColumnIndex col_i, model::ColumnIndex col_k) const {
    long double zeros_sum = 0;
    for (size_t i = 0; i < domains_[col_i]; i++) {
        zeros_sum += std::count_if(n_i_j[i].begin(), n_i_j[i].end(),
                                   [](long double val) { return val == 0; });
    }
    return zeros_sum > min_structural_zeroes_proportion_ * domains_[col_i] * domains_[col_k];
}

/* Similar to formulae (1) from "CORDS: Automatic Discovery of Correlations and Soft
   Functional Dependencies". Seems like formulae which is given in the paper is
   incorrect. */
long double Cords::CalculateChiSquared(std::vector<std::vector<long double>> const &n_i_j,
                                       std::vector<long double> const &n_i,
                                       std::vector<long double> const &n_j,
                                       model::ColumnIndex col_i, model::ColumnIndex col_k,
                                       long double sample_size) const {
    long double chi_squared = 0;
    for (size_t i = 0; i < domains_[col_i]; i++) {
        for (size_t j = 0; j < domains_[col_k]; j++) {
            if (n_i[i] * n_j[j] == 0) return 0;
            long double actual = n_i_j[i][j];
            long double expected = n_i[i] * n_j[j] / sample_size;
            chi_squared += (actual - expected) * (actual - expected) / (expected);
        }
    }
    return chi_squared;
}

void Cords::Init(model::ColumnIndex columns) {
    is_skewed_.resize(columns, false);
    domains_.resize(columns, 0);
}

bool Cords::ChiSquaredTest(model::ColumnIndex col_i, model::ColumnIndex col_k, Sample &smp,
                           std::vector<std::vector<long double>> const &n_i_j,
                           std::vector<long double> const &n_i,
                           std::vector<long double> const &n_j) const {
    long double chi_squared =
            CalculateChiSquared(n_i_j, n_i, n_j, col_i, col_k, smp.GetRowIndices().size());

    long double v = (domains_[col_i] - 1) * (domains_[col_k] - 1);

    boost::math::chi_squared dist(v);

    long double t = quantile(dist, 1 - max_false_positive_probability_);
    return chi_squared > t;
}

void Cords::RegisterCorrelation(model::ColumnIndex lhs_ind, model::ColumnIndex rhs_ind) {
    Column lhs_col(typed_relation_->GetSchema(),
                   typed_relation_->GetSchema()->GetColumn(lhs_ind)->GetName(), lhs_ind);
    Column rhs_col(typed_relation_->GetSchema(),
                   typed_relation_->GetSchema()->GetColumn(rhs_ind)->GetName(), rhs_ind);
    Correlation correlation_to_register(lhs_col, rhs_col);
    correlations_collection_.Register(std::move(correlation_to_register));
}

void Cords::CheckAndRegisterCorrelation(model::ColumnIndex col_i, model::ColumnIndex col_k,
                                        std::vector<model::TypedColumnData> const &data,
                                        Sample &smp, FrequencyHandler const &handler) {
    SkewHandling(col_i, col_k, data, smp, handler);

    std::vector<std::vector<long double>> n_i_j(domains_[col_i],
                                                std::vector<long double>(domains_[col_k], 0));
    std::vector<long double> n_i(domains_[col_i], 0);
    std::vector<long double> n_j(domains_[col_k], 0);

    FillContingencyTable(smp, data, col_i, col_k, n_i_j, n_i, n_j, handler);

    if (TooMuchStructuralZeros(n_i_j, col_i, col_k) ||
        ChiSquaredTest(col_i, col_k, smp, n_i_j, n_i, n_j)) {
        RegisterCorrelation(col_i, col_k);
        return;
    }
}

bool Cords::IsTrivialCase(model::ColumnIndex col_ind, FrequencyHandler const &handler,
                          size_t row_count) {
    if (handler.GetColumnCardinality(col_ind) >= (1 - minimum_cardinality_) * row_count) {
        Column target(typed_relation_->GetSchema(),
                      typed_relation_->GetSchema()->GetColumn(col_ind)->GetName(), col_ind);
        soft_keys_.push_back(std::move(target));
        return true;
    }
    if (handler.GetColumnCardinality(col_ind) == 1) {
        Column target(typed_relation_->GetSchema(),
                      typed_relation_->GetSchema()->GetColumn(col_ind)->GetName(), col_ind);
        trivial_columns_.push_back(target);
        return true;
    }
    return false;
}

unsigned long long Cords::ExecuteInternal() {
    std::vector<model::TypedColumnData> const &data = typed_relation_->GetColumnData();

    size_t row_count = data.front().GetNumRows();
    model::ColumnIndex column_count = data.size();

    Init(column_count);

    auto start_time = std::chrono::high_resolution_clock::now();
    FrequencyHandler handler(data, column_count, max_amount_of_categories_);

    SetProgress(kTotalProgressPercent);
    ToNextProgressPhase();

    std::vector<bool> is_soft_or_trivial;
    is_soft_or_trivial.reserve(column_count);
    for (model::ColumnIndex col_ind = 0; col_ind != column_count; ++col_ind)
        is_soft_or_trivial.push_back(IsTrivialCase(col_ind, handler, row_count));

    auto sort_indices_by_cardinality = [&handler](model::ColumnIndex ind1, model::ColumnIndex ind2)
            -> std::pair<model::ColumnIndex, model::ColumnIndex> {
        if (handler.GetColumnCardinality(ind2) > handler.GetColumnCardinality(ind1)) {
            return {ind2, ind1};
        }
        return {ind1, ind2};
    };

    for (model::ColumnIndex ind1 = 0; ind1 < column_count - 1; ind1++) {
        if (is_soft_or_trivial[ind1]) continue;

        for (model::ColumnIndex ind2 = ind1 + 1; ind2 < column_count; ind2++) {
            if (is_soft_or_trivial[ind2]) continue;

            auto [col_i, col_k] = sort_indices_by_cardinality(ind1, ind2);

            unsigned long long sample_size = Sample::CalculateSampleSize(
                    handler.GetColumnCardinality(col_i), handler.GetColumnCardinality(col_k),
                    max_false_positive_probability_, delta_);

            Sample smp(sample_size, row_count, col_i, col_k, data, typed_relation_->GetSchema());

            if (DetectAndRegisterSFD(smp) || only_sfd_) continue;

            CheckAndRegisterCorrelation(col_i, col_k, data, smp, handler);
        }
    }

    SetProgress(kTotalProgressPercent);
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    return elapsed_time.count();
}
}  // namespace algos
