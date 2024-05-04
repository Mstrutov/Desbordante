#include "frequency_handler.h"

#include <algorithm>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "model/table/column_index.h"
#include "model/table/typed_column_data.h"

namespace algos {
FrequencyHandler::FrequencyHandler(std::vector<model::TypedColumnData> const &data,
                                   model::ColumnIndex columns, size_t max_amount_of_categories) {
    cardinality_.resize(columns, 0);
    frequency_maps_.resize(columns);
    freq_sums_.resize(columns, 0);
    for (model::ColumnIndex col_ind = 0; col_ind < data.size(); col_ind++) {
        std::unordered_map<std::string, size_t> value_frequency_counter;
        auto const &col_data = data[col_ind];
        for (size_t row_ind = 0; row_ind < col_data.GetNumRows(); row_ind++) {
            auto &freq = value_frequency_counter[col_data.GetDataAsString(row_ind)];
            if (freq == 0) cardinality_[col_ind]++;
            ++freq;
        }
        std::vector<std::pair<std::string, size_t>> values_ordered_by_frequencies(
                value_frequency_counter.begin(), value_frequency_counter.end());

        auto cmp = [](std::pair<std::string, size_t> const &left,
                      std::pair<std::string, size_t> const &right) {
            return left.second > right.second;
        };

        std::sort(values_ordered_by_frequencies.begin(), values_ordered_by_frequencies.end(), cmp);

        for (size_t ordinal_number = 0;
             ordinal_number <
             std::min(max_amount_of_categories, values_ordered_by_frequencies.size());
             ordinal_number++) {
            auto const &[value, freq] = values_ordered_by_frequencies[ordinal_number];
            frequency_maps_[col_ind][value] = ordinal_number;
            freq_sums_[col_ind] += freq;
        }
    }
}

}  // namespace algos
