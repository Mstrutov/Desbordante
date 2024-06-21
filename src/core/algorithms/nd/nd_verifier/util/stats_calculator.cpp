#include "algorithms/nd/nd_verifier/util/stats_calculator.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <easylogging++.h>

#include "algorithms/nd/nd.h"

namespace algos::nd_verifier::util {

template <typename It>
std::shared_ptr<std::vector<size_t>> StatsCalculator::CalculateFrequencies(size_t codes_number,
                                                                           It begin, It end) {
    auto result = std::make_shared<std::vector<size_t>>();
    for (size_t code{0}; code < codes_number; ++code) {
        result->push_back(std::count(begin, end, code));
    }

    return result;
}

[[nodiscard]] std::unordered_map<std::string, size_t> StatsCalculator::GetLhsFrequencies() const {
    if (lhs_values_->size() != lhs_frequencies_->size()) {
        throw std::runtime_error(
                "(StatsCalulator): lhs_frequencies contains an incorrect number of values");
    }

    std::unordered_map<std::string, size_t> result;
    for (size_t i{0}; i < lhs_values_->size(); ++i) {
        result.emplace((*lhs_values_)[i], (*lhs_frequencies_)[i]);
    }

    return result;
}

[[nodiscard]] std::unordered_map<std::string, size_t> StatsCalculator::GetRhsFrequencies() const {
    if (rhs_values_->size() != rhs_frequencies_->size()) {
        throw std::runtime_error(
                "(StatsCalulator): rhs_frequencies contains an incorrect number of values");
    }

    std::unordered_map<std::string, size_t> result;
    for (size_t i{0}; i < rhs_values_->size(); ++i) {
        result.emplace((*rhs_values_)[i], (*rhs_frequencies_)[i]);
    }

    return result;
}

void StatsCalculator::CalculateStats() {
    lhs_frequencies_ =
            CalculateFrequencies(lhs_values_->size(), encoded_lhs_->begin(), encoded_lhs_->end());
    rhs_frequencies_ =
            CalculateFrequencies(rhs_values_->size(), encoded_rhs_->begin(), encoded_rhs_->end());

    std::vector<size_t> highlights_lhs_codes;
    model::WeightType max_weight{0};
    model::WeightType min_weight{UINT_MAX};
    for (size_t code{0}; code < lhs_values_->size(); ++code) {
        size_t weight = value_deps_[code].size();

        if (weight == 0) {
            // This value is associated with no rhs values, but is somehow present in value_deps
            LOG(INFO) << "WARNING: Lhs value '" << (*lhs_values_)[code]
                      << "''s weight is 0, but it's present in value_deps";
        }

        if (weight > max_weight) {
            max_weight = weight;
            highlights_lhs_codes.clear();
            highlights_lhs_codes.push_back(code);
        } else if (weight == max_weight) {
            highlights_lhs_codes.push_back(code);
        }

        if (weight < min_weight) {
            min_weight = weight;
        }
    }

    for (auto highlight_lhs_code : highlights_lhs_codes) {
        std::unordered_set<size_t>& rhs_set = value_deps_[highlight_lhs_code];

        highlights_.emplace_back(lhs_values_, rhs_values_, encoded_lhs_, encoded_rhs_,
                                 lhs_frequencies_, rhs_frequencies_, highlight_lhs_code,
                                 std::move(rhs_set));
    }

    real_weight_ = max_weight;
    global_min_weight_ = min_weight;

    LOG(INFO) << "Minimal weight: " << global_min_weight_;
    LOG(INFO) << "Maximal weight (real weight): " << real_weight_;
}

}  // namespace algos::nd_verifier::util
