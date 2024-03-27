#include "sfd_algorithm.h"

#include <chrono>
#include <string>
#include <utility>
#include <vector>

#include <boost/math/distributions/chi_squared.hpp>

#include "config/names_and_descriptions.h"
#include "config/option.h"
#include "config/tabular_data/input_table/option.h"
#include "frequency_handler.h"
#include "index_pair.h"
#include "model/table/column_index.h"
#include "model/table/typed_column_data.h"
#include "sample.h"

namespace algos {
SFDAlgorithm::SFDAlgorithm() : Algorithm({}) {
    RegisterOptions();
    MakeOptionsAvailable({config::kTableOpt.GetName()});
}

void SFDAlgorithm::RegisterOptions() {
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
        if (param <= 0) throw config::ConfigurationError("Parameter out of range");
    };
    auto check_delta = [this](long double param) {
        if (param < 0 || param > 1) throw config::ConfigurationError("delta out of range");
        if (param < minimum_cardinality_) {
            throw config::ConfigurationError("delta must be less than minimum_cardinality_");
        }
    };
    RegisterOption(config::kTableOpt(&input_table_));
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
    RegisterOption(Option{&min_structural_zeroes_amount_, kMinStructuralZeroesAmount,
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

void SFDAlgorithm::MakeExecuteOptsAvailable() {
    using namespace config::names;
    MakeOptionsAvailable({kOnlySFD, kMinCard, kMaxDiffValsProportion, kMinSFDStrengthMeasure,
                          kMinSkewThreshold, kMinStructuralZeroesAmount,
                          kMaxFalsePositiveProbability, kDelta, kMaxAmountOfCategories});
}

void SFDAlgorithm::ResetState() {
    is_skewed_.clear();
    domains_.clear();
    soft_keys_.clear();
    trivial_columns_.clear();
    correlations_.clear();
    unrelated_.clear();
    sfd_.clear();
}

void SFDAlgorithm::LoadDataInternal() {
    typed_relation_ = model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_, false);
}

bool SFDAlgorithm::DetectAndEmplaceSfd(Sample const &smp) {
    if (smp.GetConcatCardinality() <= max_diff_vals_proportion_ * smp.GetRowIndices().size()) {
        if (smp.GetLhsCardinality() >=
            (1 - min_sfd_strength_measure_) * smp.GetConcatCardinality()) {
            sfd_.emplace_back(smp.GetLhsInd(), smp.GetRhsInd());
            return true;
        } else if (smp.GetRhsCardinality() >=
                   (1 - min_sfd_strength_measure_) * smp.GetConcatCardinality()) {
            sfd_.emplace_back(smp.GetRhsInd(), smp.GetLhsInd());
            return true;
        }
    }
    return false;
}

void SFDAlgorithm::SkewHandling(model::ColumnIndex col_i, model::ColumnIndex col_k,
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

size_t SFDAlgorithm::Category(model::ColumnIndex col_ind, std::string const &val, size_t domain,
                              bool skew, FrequencyHandler const &handler) const {
    if (skew) {
        return handler.GetValueOrdinalNumberAtColumn(val, col_ind);
    }
    return std::hash<std::string>{}(val) % domain;  //+1;
}

void SFDAlgorithm::FillContingencyTable(Sample &smp,
                                        std::vector<model::TypedColumnData> const &data,
                                        model::ColumnIndex col_i, model::ColumnIndex col_k,
                                        std::vector<std::vector<long double>> &n_i_j,
                                        std::vector<long double> &n_i,
                                        std::vector<long double> &n_j,
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

bool SFDAlgorithm::TooMuchStructuralZeros(std::vector<std::vector<long double>> const &n_i_j,
                                          model::ColumnIndex col_i,
                                          model::ColumnIndex col_k) const {
    long double zeros_sum = 0;
    for (size_t i = 0; i < domains_[col_i]; i++) {
        zeros_sum += std::count_if(n_i_j[i].begin(), n_i_j[i].end(),
                                   [](long double val) { return val == 0; });
    }
    return zeros_sum > min_structural_zeroes_amount_ * domains_[col_i] * domains_[col_k];
}

/* Similar to formulae (1) from "CORDS: Automatic Discovery of Correlations and Soft
   Functional Dependencies". Seems like formulae which is given in the paper is
   incorrect. */
long double SFDAlgorithm::CalculateChiSquared(std::vector<std::vector<long double>> const &n_i_j,
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

void SFDAlgorithm::Init(model::ColumnIndex columns) {
    is_skewed_.resize(columns, false);
    domains_.resize(columns, 0);
}

std::pair<model::ColumnIndex, model::ColumnIndex> SFDAlgorithm::SwitchIndicesIfNecessary(
        model::ColumnIndex ind1, model::ColumnIndex ind2, size_t card1, size_t card2) {
    if (card2 > card1) {
        return {ind2, ind1};
    }
    return {ind1, ind2};
}

bool SFDAlgorithm::ChiSquaredTest(model::ColumnIndex col_i, model::ColumnIndex col_k, Sample &smp,
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

void SFDAlgorithm::CheckForCorrelation(model::ColumnIndex col_i, model::ColumnIndex col_k,
                                       std::vector<model::TypedColumnData> const &data, Sample &smp,
                                       FrequencyHandler const &handler) {
    SkewHandling(col_i, col_k, data, smp, handler);

    std::vector<std::vector<long double>> n_i_j(domains_[col_i],
                                                std::vector<long double>(domains_[col_k], 0));
    std::vector<long double> n_i(domains_[col_i], 0);
    std::vector<long double> n_j(domains_[col_k], 0);

    FillContingencyTable(smp, data, col_i, col_k, n_i_j, n_i, n_j, handler);

    if (TooMuchStructuralZeros(n_i_j, col_i, col_k) ||
        ChiSquaredTest(col_i, col_k, smp, n_i_j, n_i, n_j)) {
        correlations_.emplace_back(col_i, col_k);
        return;
    }
    unrelated_.emplace_back(col_i, col_k);
}

bool SFDAlgorithm::IsTrivialCase(model::ColumnIndex col_ind, std::vector<bool> &is_soft_or_trivial,
                                 FrequencyHandler const &handler, size_t row_count) {
    if (is_soft_or_trivial[col_ind]) {
        return true;
    }
    if (handler.GetColumnCardinality(col_ind) >= (1 - minimum_cardinality_) * row_count) {
        soft_keys_.push_back(col_ind);
        is_soft_or_trivial[col_ind] = true;
        return true;
    }
    if (handler.GetColumnCardinality(col_ind) == 1) {
        trivial_columns_.push_back(col_ind);
        is_soft_or_trivial[col_ind] = true;
        return true;
    }
    return false;
}

unsigned long long SFDAlgorithm::ExecuteInternal() {
    std::vector<model::TypedColumnData> const &data = typed_relation_->GetColumnData();

    size_t row_count = data.front().GetNumRows();
    model::ColumnIndex column_count = data.size();
    std::vector<bool> is_soft_or_trivial(column_count, false);
    Init(column_count);

    auto start_time = std::chrono::high_resolution_clock::now();
    FrequencyHandler handler(data, column_count, max_amount_of_categories_);

    for (model::ColumnIndex ind1 = 0; ind1 < column_count - 1; ind1++) {
        if (IsTrivialCase(ind1, is_soft_or_trivial, handler, row_count)) continue;

        for (model::ColumnIndex ind2 = ind1 + 1; ind2 < column_count; ind2++) {
            if (IsTrivialCase(ind2, is_soft_or_trivial, handler, row_count)) continue;

            auto [col_i, col_k] =
                    SwitchIndicesIfNecessary(ind1, ind2, handler.GetColumnCardinality(ind1),
                                             handler.GetColumnCardinality(ind2));

            unsigned long long sample_size = Sample::CalculateSampleSize(
                    handler.GetColumnCardinality(col_i), handler.GetColumnCardinality(col_k),
                    max_false_positive_probability_, delta_);

            Sample smp(sample_size, row_count, col_i, col_k, data);

            if (DetectAndEmplaceSfd(smp) || only_sfd_) continue;

            CheckForCorrelation(col_i, col_k, data, smp, handler);
        }
    }
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);
    return elapsed_time.count();
}

std::vector<model::ColumnIndex> const &SFDAlgorithm::GetSoftKeys() const {
    return soft_keys_;
}

std::vector<model::ColumnIndex> const &SFDAlgorithm::GetTrivialColumns() const {
    return trivial_columns_;
}

std::vector<ColumnIndexPair> const &SFDAlgorithm::GetCorrelations() const {
    return correlations_;
}

std::vector<ColumnIndexPair> const &SFDAlgorithm::GetUnrelated() const {
    return unrelated_;
}

std::vector<ColumnIndexPair> const &SFDAlgorithm::GetSFD() const {
    return sfd_;
}
}  // namespace algos
