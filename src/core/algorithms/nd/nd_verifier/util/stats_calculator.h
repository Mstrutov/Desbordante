#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "algorithms/nd/nd.h"
#include "algorithms/nd/nd_verifier/util/highlight.h"

namespace algos::nd_verifier::util {

class StatsCalculator {
private:
    std::unordered_map<size_t, std::unordered_set<size_t>> value_deps_;

    // Shared data:
    std::shared_ptr<std::vector<std::string>> lhs_values_;
    std::shared_ptr<std::vector<std::string>> rhs_values_;
    std::shared_ptr<std::vector<size_t>> encoded_lhs_;
    std::shared_ptr<std::vector<size_t>> encoded_rhs_;

    // Cached data:
    std::vector<Highlight> highlights_;
    std::shared_ptr<std::vector<size_t>> lhs_frequencies_;
    std::shared_ptr<std::vector<size_t>> rhs_frequencies_;
    model::WeightType global_min_weight_{UINT_MAX};
    model::WeightType real_weight_{0};

    template <typename It>
    std::shared_ptr<std::vector<size_t>> CalculateFrequencies(size_t codes_number, It begin,
                                                              It end);

public:
    StatsCalculator(std::unordered_map<size_t, std::unordered_set<size_t>>&& value_deps,
                    std::shared_ptr<std::vector<std::string>> lhs_codes,
                    std::shared_ptr<std::vector<std::string>> rhs_codes,
                    std::shared_ptr<std::vector<size_t>> encoded_lhs,
                    std::shared_ptr<std::vector<size_t>> encoded_rhs)
        : value_deps_(std::move(value_deps)),
          lhs_values_(std::move(lhs_codes)),
          rhs_values_(std::move(rhs_codes)),
          encoded_lhs_(std::move(encoded_lhs)),
          encoded_rhs_(std::move(encoded_rhs)) {}

    StatsCalculator() = default;

    [[nodiscard]] std::vector<Highlight> const& GetHighlights() const {
        return highlights_;
    }

    [[nodiscard]] model::WeightType GetGlobalMinWeight() const {
        return global_min_weight_;
    }

    [[nodiscard]] model::WeightType GetRealWeight() const {
        return real_weight_;
    }

    [[nodiscard]] std::unordered_map<std::string, size_t> GetLhsFrequencies() const;

    [[nodiscard]] std::unordered_map<std::string, size_t> GetRhsFrequencies() const;

    void CalculateStats();
};

}  // namespace algos::nd_verifier::util
