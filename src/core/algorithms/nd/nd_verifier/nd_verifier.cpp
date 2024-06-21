#include "algorithms/nd/nd_verifier/nd_verifier.h"

#include <chrono>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <easylogging++.h>

#include "algorithms/nd/nd_verifier/util/stats_calculator.h"
#include "algorithms/nd/nd_verifier/util/vector_to_string.h"
#include "config/descriptions.h"
#include "config/equal_nulls/option.h"
#include "config/indices/option.h"
#include "config/names.h"
#include "config/option.h"
#include "config/option_using.h"
#include "config/tabular_data/input_table/option.h"
#include "model/table/column_layout_typed_relation_data.h"
#include "model/table/typed_column_data.h"
#include "model/types/type.h"
#include "util/timed_invoke.h"

namespace algos::nd_verifier {

NDVerifier::NDVerifier() : Algorithm({}) {
    RegisterOptions();
    MakeOptionsAvailable({config::kTableOpt.GetName(), config::kEqualNullsOpt.GetName()});
}

void NDVerifier::RegisterOptions() {
    DESBORDANTE_OPTION_USING;

    auto get_schema_cols = [this]() { return typed_relation_->GetSchema()->GetNumColumns(); };

    RegisterOption(config::kTableOpt(&input_table_));
    RegisterOption(config::kEqualNullsOpt(&is_null_equal_null_));
    RegisterOption(config::kLhsIndicesOpt(&lhs_indices_, get_schema_cols));
    RegisterOption(config::kRhsIndicesOpt(&rhs_indices_, get_schema_cols));
    RegisterOption(Option<model::WeightType>(&weight_, kWeight, kDNDWeight, 1));
}

void NDVerifier::LoadDataInternal() {
    typed_relation_ =
            model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_, is_null_equal_null_);
    input_table_->Reset();
    if (typed_relation_->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty dataset: ND verifying is meaningless.");
    }
}

void NDVerifier::MakeExecuteOptsAvailable() {
    MakeOptionsAvailable({config::kLhsIndicesOpt.GetName(), config::kRhsIndicesOpt.GetName(),
                          config::names::kWeight});
}

unsigned long long NDVerifier::ExecuteInternal() {
    LOG(INFO) << "Parameters of NDVerifier:";
    LOG(INFO) << "\tInput table: " << input_table_->GetRelationName();
    LOG(INFO) << "\tNull equals null: " << is_null_equal_null_;
    LOG(INFO) << "\tLhs indices: " << util::VectorToString(lhs_indices_);
    LOG(INFO) << "\tRhs indices: " << util::VectorToString(rhs_indices_);
    LOG(INFO) << "\tWeight: " << weight_;

    auto verification_time = ::util::TimedInvoke(&NDVerifier::VerifyND, this);

    LOG(DEBUG) << "ND verification took " << std::to_string(verification_time) << "ms";

    auto stats_calculation_time = ::util::TimedInvoke(&NDVerifier::CalculateStats, this);

    LOG(DEBUG) << "Statistics calculation took " << std::to_string(stats_calculation_time) << "ms";

    return verification_time + stats_calculation_time;
}

bool NDVerifier::NDHolds() const {
    return stats_calculator_.GetRealWeight() <= weight_;
}

void NDVerifier::ResetState() {
    stats_calculator_ = util::StatsCalculator{};
}

NDVerifier::EncodedValuesType NDVerifier::EncodeMultipleValues(
        config::IndicesType const& col_idxs) const {
    auto values = std::make_shared<std::vector<std::string>>();
    auto row = std::make_shared<std::vector<size_t>>();

    for (size_t row_idx{0}; row_idx < typed_relation_->GetNumRows(); ++row_idx) {
        std::stringstream ss;
        ss << '(';
        bool is_null = false;
        for (auto col_idx_pt{col_idxs.begin()}; col_idx_pt != col_idxs.end(); ++col_idx_pt) {
            model::TypedColumnData const& col_data = typed_relation_->GetColumnData(*col_idx_pt);
            std::vector<std::byte const*> const& byte_data = col_data.GetData();
            model::Type const& type = col_data.GetType();

            if (col_idx_pt != col_idxs.begin()) {
                ss << ", ";
            }

            std::byte const* bytes_ptr = byte_data[row_idx];
            if (bytes_ptr == nullptr) {
                LOG(INFO) << "WARNING: Cell (" << *col_idx_pt << ", " << row_idx << ") is empty";
                is_null = true;
            } else {
                ss << type.ValueToString(bytes_ptr);
            }
        }
        ss << ')';

        auto string_data = ss.str();
        size_t index{values->size()};
        if (is_null && !is_null_equal_null_) {
            // If not null_eq_null, every value with null will be unique
            values->push_back(string_data);
        } else {
            for (size_t i{0}; i < values->size(); ++i) {
                if ((*values)[i] == string_data) {
                    index = i;
                    break;
                }
            }

            if (index == values->size()) {  // wasn't in codes
                values->push_back(string_data);
            }
        }
        row->push_back(index);
    }
    return std::make_pair(std::move(values), std::move(row));
}

NDVerifier::EncodedValuesType NDVerifier::EncodeSingleValues(config::IndexType col_idx) const {
    auto values = std::make_shared<std::vector<std::string>>();
    auto row = std::make_shared<std::vector<size_t>>();

    model::TypedColumnData const& col_data = typed_relation_->GetColumnData(col_idx);
    std::vector<std::byte const*> const& byte_data = col_data.GetData();
    model::Type const& type = col_data.GetType();

    for (std::byte const* bytes_ptr : byte_data) {
        std::string string_data;
        bool is_null = false;
        if (bytes_ptr == nullptr) {
            LOG(INFO) << "WARNING: Empty cell in column " << col_idx;
            is_null = true;
        } else {
            string_data = type.ValueToString(bytes_ptr);
        }

        size_t index{values->size()};
        if (is_null && !is_null_equal_null_) {
            // If not null_eq_null, every value with null will be unique
            values->push_back(string_data);
        } else {
            for (size_t i{0}; i < values->size(); ++i) {
                if ((*values)[i] == string_data) {
                    index = i;
                    break;
                }
            }

            if (index == values->size()) {  // wasn`t in values
                values->push_back(string_data);
            }
        }
        row->push_back(index);
    }
    return std::make_pair(std::move(values), std::move(row));
}

NDVerifier::EncodedValuesType NDVerifier::EncodeValues(config::IndicesType const& col_idxs) const {
    return col_idxs.size() == 1 ? EncodeSingleValues(col_idxs[0]) : EncodeMultipleValues(col_idxs);
}

void NDVerifier::VerifyND() {
    auto local_start_time = std::chrono::system_clock::now();
    auto [lhs_values, encoded_lhs] = EncodeValues(lhs_indices_);
    auto [rhs_values, encoded_rhs] = EncodeValues(rhs_indices_);

    LOG(DEBUG) << "Values encoding took "
               << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                                         std::chrono::system_clock::now() - local_start_time)
                                         .count())
               << "ms";

    local_start_time = std::chrono::system_clock::now();
    std::unordered_map<size_t, std::unordered_set<size_t>> value_deps{};

    for (size_t i{0}; i < encoded_lhs->size(); ++i) {
        auto lhs_code = (*encoded_lhs)[i];
        auto rhs_code = (*encoded_rhs)[i];

        if (value_deps.find(lhs_code) == value_deps.end()) {
            value_deps.emplace(lhs_code, std::unordered_set<size_t>({rhs_code}));
        } else {
            value_deps[lhs_code].insert(rhs_code);
        }
    }
    LOG(DEBUG) << "Value deps calculation took "
               << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                                         std::chrono::system_clock::now() - local_start_time)
                                         .count())
               << "ms";

    stats_calculator_ = util::StatsCalculator(std::move(value_deps), std::move(lhs_values),
                                              std::move(rhs_values), std::move(encoded_lhs),
                                              std::move(encoded_rhs));
}

[[nodiscard]] std::vector<util::Highlight> const& NDVerifier::GetHighlights() const {
    return stats_calculator_.GetHighlights();
}

[[nodiscard]] model::WeightType NDVerifier::GetGlobalMinWeight() const {
    return stats_calculator_.GetGlobalMinWeight();
}

[[nodiscard]] model::WeightType NDVerifier::GetRealWeight() const {
    return stats_calculator_.GetRealWeight();
}

[[nodiscard]] std::unordered_map<std::string, size_t> NDVerifier::GetLhsFrequencies() const {
    return stats_calculator_.GetLhsFrequencies();
}

[[nodiscard]] std::unordered_map<std::string, size_t> NDVerifier::GetRhsFrequencies() const {
    return stats_calculator_.GetRhsFrequencies();
}

}  // namespace algos::nd_verifier
