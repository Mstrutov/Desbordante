#pragma once
#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "frequency_handler.h"
#include "model/table/column_index.h"
#include "model/table/typed_column_data.h"

namespace algos {

class Sample {
private:
    std::vector<size_t> row_indices_;
    model::ColumnIndex lhs_ind_;
    model::ColumnIndex rhs_ind_;
    size_t lhs_cardinality_;
    size_t rhs_cardinality_;
    size_t concat_cardinality_;

public:
    Sample(unsigned long long sample_size, size_t rows, model::ColumnIndex l, model::ColumnIndex r,
           std::vector<model::TypedColumnData> const &data);
    void Filter(FrequencyHandler const &handler, std::vector<model::TypedColumnData> const &data,
                model::ColumnIndex col_ind);

    /* Formulae (2) from "CORDS: Automatic Discovery of Correlations and Soft Functional
       Dependencies."*/
    static long long CalculateSampleSize(size_t d1, size_t d2,
                                         long double max_false_positive_probability,
                                         long double delta) {
        long double v = (d1 - 1) * (d2 - 1);
        long double d = std::min(d1, d2);
        long double log =
                std::log(max_false_positive_probability * std::sqrt(2 * std::numbers::pi));
        long double numerator = std::pow(-16 * v * log, 0.5) - 8 * log;
        long double denominator = delta * (d - 1);
        long double v2 = std::pow(v, 0.071);
        return static_cast<long long>((numerator / denominator) * (v2 / 1.69));
    }

    [[nodiscard]] std::vector<size_t> const &GetRowIndices() const;
    [[nodiscard]] model::ColumnIndex GetLhsInd() const;
    [[nodiscard]] model::ColumnIndex GetRhsInd() const;
    [[nodiscard]] model::ColumnIndex GetLhsCardinality() const;
    [[nodiscard]] model::ColumnIndex GetRhsCardinality() const;
    [[nodiscard]] model::ColumnIndex GetConcatCardinality() const;
};
}  // namespace algos
