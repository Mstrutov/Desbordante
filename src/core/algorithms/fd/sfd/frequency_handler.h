#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "model/table/column_index.h"
#include "model/table/typed_column_data.h"

namespace algos {

class FrequencyHandler {
private:
    std::vector<size_t> cardinality_;
    std::vector<size_t> freq_sums_;
    std::vector<std::unordered_map<std::string, size_t>> frequency_maps_;

public:
    FrequencyHandler(std::vector<model::TypedColumnData> const &data, model::ColumnIndex columns,
                     size_t max_amount_of_categories);

    [[nodiscard]] size_t GetColumnFrequencySum(model::ColumnIndex col_ind) const {
        return freq_sums_[col_ind];
    }

    [[nodiscard]] size_t GetColumnCardinality(model::ColumnIndex col_ind) const {
        return cardinality_[col_ind];
    }

    [[nodiscard]] size_t GetValueOrdinalNumberAtColumn(std::string const &val,
                                                       model::ColumnIndex col_ind) const {
        return frequency_maps_[col_ind].at(val);
    }

    [[nodiscard]] bool ContainsValAtColumn(std::string const &val,
                                           model::ColumnIndex col_ind) const {
        return frequency_maps_[col_ind].find(val) != frequency_maps_[col_ind].end();
    }

    [[nodiscard]] size_t Size() const {
        return frequency_maps_.size();
    }

    [[nodiscard]] size_t ColumnFrequencyMapSize(model::ColumnIndex col_ind) const {
        return frequency_maps_[col_ind].size();
    }

    void Clear() {
        cardinality_.clear();
        frequency_maps_.clear();
        freq_sums_.clear();
    }
};
}  // namespace algos
