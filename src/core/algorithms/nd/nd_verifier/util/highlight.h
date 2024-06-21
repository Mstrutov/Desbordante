#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <unordered_set>
#include <vector>

namespace algos::nd_verifier::util {

/// @brief Lhs value that violates ND
class Highlight {
private:
    // Shared data:
    std::shared_ptr<std::vector<std::string>> lhs_values_;
    std::shared_ptr<std::vector<std::string>> rhs_values_;
    std::shared_ptr<std::vector<size_t>> encoded_lhs_;
    std::shared_ptr<std::vector<size_t>> encoded_rhs_;
    std::shared_ptr<std::vector<size_t>> lhs_frequencies_;
    std::shared_ptr<std::vector<size_t>> rhs_frequencies_;

    // Highlight parameters:
    size_t highlight_lhs_;
    std::unordered_set<size_t> highlight_rhs_;

    // Cached data:
    std::unique_ptr<std::unordered_set<size_t>> most_frequent_rhs_value_codes_{nullptr};
    std::unique_ptr<std::vector<size_t>> occurences_indices_{nullptr};
    std::unique_ptr<std::vector<std::string>> highlight_rhs_values_{nullptr};

    /// @brief Find most frequent rhs values and cache them
    std::unordered_set<size_t> const& GetMostFrequentRhsValueCodes();

public:
    Highlight(std::shared_ptr<std::vector<std::string>> lhs_codes,
              std::shared_ptr<std::vector<std::string>> rhs_codes,
              std::shared_ptr<std::vector<size_t>> encoded_lhs,
              std::shared_ptr<std::vector<size_t>> encoded_rhs,
              std::shared_ptr<std::vector<size_t>> lhs_frequencies,
              std::shared_ptr<std::vector<size_t>> rhs_frequencies, size_t highlight_lhs_code,
              std::unordered_set<size_t>&& highlight_rhs_codes)
        : lhs_values_(std::move(lhs_codes)),
          rhs_values_(std::move(rhs_codes)),
          encoded_lhs_(std::move(encoded_lhs)),
          encoded_rhs_(std::move(encoded_rhs)),
          lhs_frequencies_(std::move(lhs_frequencies)),
          rhs_frequencies_(std::move(rhs_frequencies)),
          highlight_lhs_(highlight_lhs_code),
          highlight_rhs_(std::move(highlight_rhs_codes)) {}

    /// @brief Get cached indices of all occurences of this highlight of calculate it (heavy
    /// operation)
    [[nodiscard]] std::vector<size_t> const& GetOccurencesIndices();

    /// @brief Get number of occurences of this highlight
    [[nodiscard]] size_t GetOccurencesNumber();

    /// @brief Get Lhs value
    [[nodiscard]] std::string GetLhsValue() const {
        return (*lhs_values_)[highlight_lhs_];
    }

    /// @brief Get cached Rhs values or calculate it (heavy operation)
    [[nodiscard]] std::vector<std::string> const& GetRhsValues();

    /// @brief Get number of distinct values in rhs
    [[nodiscard]] size_t GetDistinctRhsValuesNumber() const {
        return highlight_rhs_.size();
    }

    /// @brief Get indices of most frequent rhs values -- heavy operation
    [[nodiscard]] std::vector<size_t> GetMostFrequentRhsValueIndices();

    /// @brief Get most frequent rhs values -- heavy operation
    [[nodiscard]] std::vector<std::string> GetMostFrequentRhsValues();

    [[nodiscard]] std::string ToIndicesString();

    [[nodiscard]] std::string ToValuesString();
};

std::ostream& operator<<(std::ostream& os, Highlight& hl);

}  // namespace algos::nd_verifier::util
