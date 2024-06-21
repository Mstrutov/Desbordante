#include "algorithms/nd/nd_verifier/util/highlight.h"

#include <algorithm>
#include <memory>
#include <ostream>
#include <string>
#include <unordered_set>
#include <vector>

#include "algorithms/nd/nd_verifier/util/vector_to_string.h"

namespace algos::nd_verifier::util {

std::unordered_set<size_t> const& Highlight::GetMostFrequentRhsValueCodes() {
    if (most_frequent_rhs_value_codes_ == nullptr) {
        most_frequent_rhs_value_codes_ = std::make_unique<std::unordered_set<size_t>>();
        size_t num_most_frequent_rhs_values{0};
        for (size_t rhs_code : highlight_rhs_) {
            size_t freq = (*rhs_frequencies_)[rhs_code];
            if (freq > num_most_frequent_rhs_values) {
                num_most_frequent_rhs_values = freq;
                most_frequent_rhs_value_codes_->clear();
                most_frequent_rhs_value_codes_->insert(rhs_code);
            } else if (freq == num_most_frequent_rhs_values) {
                most_frequent_rhs_value_codes_->insert(rhs_code);
            }
        }
    }

    return *most_frequent_rhs_value_codes_;
}

[[nodiscard]] std::vector<size_t> const& Highlight::GetOccurencesIndices() {
    if (occurences_indices_ == nullptr) {
        occurences_indices_ = std::make_unique<std::vector<size_t>>();
        for (size_t i{0}; i < encoded_lhs_->size(); ++i) {
            if ((*encoded_lhs_)[i] == highlight_lhs_) {
                occurences_indices_->push_back(i);
            }
        }
    }

    return *occurences_indices_;
}

[[nodiscard]] size_t Highlight::GetOccurencesNumber() {
    return GetOccurencesIndices().size();
}

[[nodiscard]] std::vector<std::string> const& Highlight::GetRhsValues() {
    if (highlight_rhs_values_ == nullptr) {
        highlight_rhs_values_ = std::make_unique<std::vector<std::string>>();
        std::transform(highlight_rhs_.begin(), highlight_rhs_.end(),
                       std::back_inserter(*highlight_rhs_values_),
                       [this](size_t code) { return (*rhs_values_)[code]; });
    }

    return *highlight_rhs_values_;
}

[[nodiscard]] std::vector<size_t> Highlight::GetMostFrequentRhsValueIndices() {
    std::vector<size_t> most_freq_rhs_inds;

    std::unordered_set<size_t> const& most_freq_rhs_codes = GetMostFrequentRhsValueCodes();
    for (size_t idx{0}; idx < encoded_rhs_->size(); ++idx) {
        size_t candidate_code = (*encoded_rhs_)[idx];
        if (most_freq_rhs_codes.find(candidate_code) != most_freq_rhs_codes.end()) {
            most_freq_rhs_inds.push_back(idx);
        }
    }

    return most_freq_rhs_inds;
}

[[nodiscard]] std::vector<std::string> Highlight::GetMostFrequentRhsValues() {
    std::vector<std::string> most_freq_rhs_vals;
    std::unordered_set<size_t> const& most_freq_rhs_codes = GetMostFrequentRhsValueCodes();

    std::transform(most_freq_rhs_codes.begin(), most_freq_rhs_codes.end(),
                   std::back_inserter(most_freq_rhs_vals),
                   [this](size_t code) { return (*rhs_values_)[code]; });

    return most_freq_rhs_vals;
}

[[nodiscard]] std::string Highlight::ToIndicesString() {
    return util::VectorToString(GetOccurencesIndices());
}

[[nodiscard]] std::string Highlight::ToValuesString() {
    std::string const& lhs = GetLhsValue();
    std::vector<std::string> const& rhs = GetRhsValues();

    return lhs + " -> " + util::VectorToString(rhs);
}

std::ostream& operator<<(std::ostream& os, Highlight& hl) {
    os << hl.ToValuesString();
    return os;
}

}  // namespace algos::nd_verifier::util
