#include "algorithms/nd/nd_verifier/nd_verifier.h"

#include <algorithm>
#include <chrono>
#include <set>
#include <sstream>

#include <easylogging++.h>

#include "config/descriptions.h"
#include "config/equal_nulls/option.h"
#include "config/indices/option.h"
#include "config/names.h"
#include "config/option.h"
#include "config/option_using.h"
#include "config/tabular_data/input_table/option.h"
#include "nd_verifier.h"

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
    typed_relation_ =
            model::ColumnLayoutTypedRelationData::CreateFrom(*input_table_, is_null_equal_null_);
}

void NDVerifier::MakeExecuteOptsAvailable() {
    MakeOptionsAvailable({config::kLhsIndicesOpt.GetName(), config::kRhsIndicesOpt.GetName(),
                          config::names::kWeight});
}

unsigned long long NDVerifier::ExecuteInternal() {
    auto start_time = std::chrono::system_clock::now();

    auto const& indices_to_str = [](config::IndicesType const& vect) -> std::string {
        return std::accumulate(std::next(vect.begin()), vect.end(), std::to_string(vect[0]),
                               [](std::string&& a, config::IndexType b) {
                                   return std::move(a) + ", " + std::to_string(b);
                               });
    };

    LOG(INFO) << "Parameters of NDVerifier:";
    LOG(INFO) << "\tInput table: " << input_table_->GetRelationName();
    LOG(INFO) << "\tNull equals null: " << is_null_equal_null_;
    LOG(INFO) << "\tLhs indices: " << indices_to_str(lhs_indices_);
    LOG(INFO) << "\tRhs indices: " << indices_to_str(rhs_indices_);
    LOG(INFO) << "\tWeight: " << weight_;

    VerifyND();

    auto elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);

    // TODO(senichenkov): here should be INFO:
    LOG(WARNING) << "NDVerifier::ExecuteInternal finished in "
                 << std::to_string(elapsed_milliseconds.count())
                 << "ms";  // We use std::to_string, because compiler on github doesn`t
                           // have implementation for stringstream::operator<<(unsigned)

    CalculateStats();

    elapsed_milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now() - start_time);

    return elapsed_milliseconds.count();
}

bool NDVerifier::NDHolds() {
    return stats_calculator_.GetRealWeight() <= weight_;
}

void NDVerifier::ResetState() {
    stats_calculator_ = util::StatsCalculator<std::string>{};
}

std::pair<std::shared_ptr<std::vector<std::string>>, std::shared_ptr<std::vector<size_t>>>
NDVerifier::EncodeMultipleValues(config::IndicesType const& col_idxs) const {
    auto values = std::make_shared<std::vector<std::string>>();
    auto row = std::make_shared<std::vector<size_t>>();

    for (size_t row_idx{0}; row_idx < typed_relation_->GetNumRows(); ++row_idx) {
        std::stringstream ss;
        ss << '(';
        bool is_null = false;
        for (auto col_idx_pt{col_idxs.begin()}; col_idx_pt != col_idxs.end(); ++col_idx_pt) {
            auto const& col_data = typed_relation_->GetColumnData(*col_idx_pt);
            auto const& byte_data = col_data.GetData();
            auto const& type = col_data.GetType();

            if (col_idx_pt != col_idxs.begin()) {
                ss << ", ";
            }

            auto const* bytes_ptr = byte_data[row_idx];
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
                if (values->at(i) == string_data) {
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

std::pair<std::shared_ptr<std::vector<std::string>>, std::shared_ptr<std::vector<size_t>>>
NDVerifier::EncodeSingleValues(config::IndexType col_idx) const {
    LOG(INFO) << "Encoding single values...";
    auto values = std::make_shared<std::vector<std::string>>();
    auto row = std::make_shared<std::vector<size_t>>();

    auto const& col_data = typed_relation_->GetColumnData(col_idx);
    auto const& byte_data = col_data.GetData();
    auto const& type = col_data.GetType();

    for (auto const* bytes_ptr : byte_data) {
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
                if (values->at(i) == string_data) {
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
    LOG(INFO) << "\t\tDone";
    return std::make_pair(std::move(values), std::move(row));
}

std::pair<std::shared_ptr<std::vector<std::string>>, std::shared_ptr<std::vector<size_t>>>
NDVerifier::EncodeValues(config::IndicesType const& col_idxs) const {
    if (col_idxs.size() == 1) {
        return EncodeSingleValues(col_idxs[0]);
    }
    return EncodeMultipleValues(col_idxs);
}

void NDVerifier::VerifyND() {
    auto local_start_time = std::chrono::system_clock::now();
    LOG(INFO) << "VerifyND started";
    auto [lhs_values, encoded_lhs] = EncodeValues(lhs_indices_);
    LOG(INFO) << "Lhs encoded";
    auto [rhs_values, encoded_rhs] = EncodeValues(rhs_indices_);
    LOG(INFO) << "Rhs encoded";

    LOG(WARNING) << "Values encoding encoding took "
                 << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::system_clock::now() - local_start_time).count());

    local_start_time = std::chrono::system_clock::now();
    auto value_deps = std::make_shared<std::unordered_map<size_t, std::unordered_set<size_t>>>();

    for (size_t i{0}; i < encoded_lhs->size(); ++i) {
        auto lhs_code = encoded_lhs->at(i);
        auto rhs_code = encoded_rhs->at(i);

        if (value_deps->find(lhs_code) == value_deps->end()) {
            value_deps->emplace(lhs_code, std::unordered_set<size_t>({rhs_code}));
        } else {
            value_deps->at(lhs_code).insert(rhs_code);
        }
    }
    LOG(WARNING) << "Value deps calculation took "
                 << std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::system_clock::now() - local_start_time).count());

    stats_calculator_ = util::StatsCalculator<std::string>(
            std::move(value_deps), std::move(lhs_values), std::move(rhs_values),
            std::move(encoded_lhs), std::move(encoded_rhs));
}

}  // namespace algos::nd_verifier
