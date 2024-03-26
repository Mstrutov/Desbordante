#include "sample.h"

#include <chrono>
#include <functional>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "frequency_handler.h"

namespace algos {
Sample::Sample(unsigned long long sample_size, size_t rows, model::ColumnIndex l,
               model::ColumnIndex r, std::vector<model::TypedColumnData> const &data)
    : lhs_ind_(l), rhs_ind_(r) {
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    std::uniform_int_distribution<model::ColumnIndex> distribution(0, rows - 1);

    std::unordered_set<std::string> map_lhs;
    std::unordered_set<std::string> map_rhs;
    std::unordered_set<std::string> map_cardinality;

    for (model::ColumnIndex i = 0; i < sample_size; i++) {
        size_t row = distribution(gen);
        row_indices_.push_back(row);
        map_lhs.insert(data[l].GetDataAsString(row));
        map_rhs.insert(data[r].GetDataAsString(row));
        map_cardinality.insert(data[l].GetDataAsString(row) + data[r].GetDataAsString(row));
    }
    lhs_cardinality_ = map_lhs.size();
    rhs_cardinality_ = map_rhs.size();
    concat_cardinality_ = map_cardinality.size();
}

void Sample::Filter(FrequencyHandler const &handler,
                    std::vector<model::TypedColumnData> const &data, model::ColumnIndex col_ind) {
    auto row = row_indices_.begin();
    while (row != row_indices_.end()) {
        if (!handler.ContainsValAtColumn(data[col_ind].GetDataAsString(*row), col_ind)) {
            row = row_indices_.erase(row);
        } else {
            row++;
        }
    }
}

[[nodiscard]] std::vector<size_t> const &Sample::GetRowIndices() const {
    return row_indices_;
}

[[nodiscard]] model::ColumnIndex Sample::GetLhsInd() const {
    return lhs_ind_;
}

[[nodiscard]] model::ColumnIndex Sample::GetRhsInd() const {
    return rhs_ind_;
}

[[nodiscard]] model::ColumnIndex Sample::GetLhsCardinality() const {
    return lhs_cardinality_;
}

[[nodiscard]] model::ColumnIndex Sample::GetRhsCardinality() const {
    return rhs_cardinality_;
}

[[nodiscard]] model::ColumnIndex Sample::GetConcatCardinality() const {
    return concat_cardinality_;
}

}  // namespace algos
