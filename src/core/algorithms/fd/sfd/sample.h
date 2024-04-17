#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>

#include "frequency_handler.h"
#include "model/table/column.h"
#include "model/table/column_index.h"
#include "model/table/relational_schema.h"
#include "model/table/typed_column_data.h"
#include "model/table/vertical.h"

namespace algos {

class Sample {
private:
    std::vector<size_t> row_indices_;
    Column lhs_col_;
    Column rhs_col_;
    size_t lhs_cardinality_;
    size_t rhs_cardinality_;
    size_t concat_cardinality_;

public:
    Sample(unsigned long long sample_size, size_t rows, model::ColumnIndex l, model::ColumnIndex r,
           std::vector<model::TypedColumnData> const &data, RelationalSchema const *rel_schema_);
    void Filter(FrequencyHandler const &handler, std::vector<model::TypedColumnData> const &data,
                model::ColumnIndex col_ind);

    /* Formulae (2) from "CORDS: Automatic Discovery of Correlations and Soft Functional
       Dependencies."*/
    static unsigned long long CalculateSampleSize(size_t d1, size_t d2,
                                                  long double max_false_positive_probability,
                                                  long double delta);

    [[nodiscard]] std::vector<size_t> const &GetRowIndices() const {
        return row_indices_;
    }

    [[nodiscard]] Vertical GetLhsVertical() const {
        return Vertical(lhs_col_);
    }

    [[nodiscard]] Column GetLhsColumn() const {
        return lhs_col_;
    }

    [[nodiscard]] Column GetRhsColumn() const {
        return rhs_col_;
    }

    [[nodiscard]] Vertical GetRhsVertical() const {
        return Vertical(rhs_col_);
    }

    [[nodiscard]] model::ColumnIndex GetLhsCardinality() const {
        return lhs_cardinality_;
    }

    [[nodiscard]] model::ColumnIndex GetRhsCardinality() const {
        return rhs_cardinality_;
    }

    [[nodiscard]] model::ColumnIndex GetConcatCardinality() const {
        return concat_cardinality_;
    }
};
}  // namespace algos
