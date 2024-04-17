#include "sample.h"

#include <chrono>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

#include "frequency_handler.h"

namespace algos {
Sample::Sample(unsigned long long sample_size, size_t rows, model::ColumnIndex l,
               model::ColumnIndex r, std::vector<model::TypedColumnData> const &data,
               RelationalSchema const *rel_schema_)
    : lhs_col_(rel_schema_, rel_schema_->GetColumn(l)->GetName(), l),
      rhs_col_(rel_schema_, rel_schema_->GetColumn(r)->GetName(), r) {
    auto seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> distribution(0, rows - 1);

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

unsigned long long Sample::CalculateSampleSize(size_t d1, size_t d2,
                                               long double max_false_positive_probability,
                                               long double delta) {
    long double v = (d1 - 1) * (d2 - 1);
    long double d = std::min(d1, d2);
    long double log = std::log(max_false_positive_probability * std::sqrt(2 * std::numbers::pi));
    long double numerator = std::pow(-16 * v * log, 0.5) - 8 * log;
    long double denominator = delta * (d - 1);
    long double v2 = std::pow(v, 0.071);
    return static_cast<long long>((numerator / denominator) * (v2 / 1.69));
}

void Sample::Filter(FrequencyHandler const &handler,
                    std::vector<model::TypedColumnData> const &data, model::ColumnIndex col_ind) {
    std::erase_if(row_indices_, [&handler, &data, col_ind](size_t row_id) {
        return !handler.ContainsValAtColumn(data[col_ind].GetDataAsString(row_id), col_ind);
    });
}
}  // namespace algos
