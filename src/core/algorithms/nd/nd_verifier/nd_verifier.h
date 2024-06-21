#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "algorithms/algorithm.h"
#include "algorithms/nd/nd.h"
#include "algorithms/nd/nd_verifier/util/highlight.h"
#include "algorithms/nd/nd_verifier/util/stats_calculator.h"
#include "config/equal_nulls/type.h"
#include "config/indices/type.h"
#include "config/tabular_data/input_table_type.h"
#include "model/table/column_layout_typed_relation_data.h"

namespace algos::nd_verifier {

/// @brief Algorithm for verifying if ND holds with given weight
class NDVerifier : public Algorithm {
private:
    using EncodedValuesType = std::pair<std::shared_ptr<std::vector<std::string>>,
                                        std::shared_ptr<std::vector<size_t>>>;

    config::InputTable input_table_;
    config::IndicesType lhs_indices_;
    config::IndicesType rhs_indices_;
    model::WeightType weight_;
    config::EqNullsType is_null_equal_null_;

    std::shared_ptr<model::ColumnLayoutTypedRelationData> typed_relation_;

    util::StatsCalculator stats_calculator_;

    void RegisterOptions();
    void ResetState() override;
    void VerifyND();

    EncodedValuesType EncodeMultipleValues(config::IndicesType const& col_idxs) const;
    EncodedValuesType EncodeSingleValues(config::IndexType col_idx) const;
    EncodedValuesType EncodeValues(config::IndicesType const& col_idxs) const;

    void CalculateStats() {
        stats_calculator_.CalculateStats();
    }

protected:
    void LoadDataInternal() override;
    void MakeExecuteOptsAvailable() override;
    unsigned long long ExecuteInternal() override;

public:
    NDVerifier();
    [[nodiscard]] bool NDHolds() const;
    [[nodiscard]] std::vector<util::Highlight> const& GetHighlights() const;
    [[nodiscard]] model::WeightType GetGlobalMinWeight() const;
    [[nodiscard]] model::WeightType GetRealWeight() const;
    [[nodiscard]] std::unordered_map<std::string, size_t> GetLhsFrequencies() const;
    [[nodiscard]] std::unordered_map<std::string, size_t> GetRhsFrequencies() const;
};

}  // namespace algos::nd_verifier
