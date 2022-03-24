//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <cmath>
#include <vector>

#include "ColumnData.h"
#include "CSVParser.h"
#include "RelationalSchema.h"
#include "RelationData.h"

class ColumnLayoutRelationData final : public RelationData {
public:
    static constexpr int kNullValueId = -1;

    using RelationData::AbstractRelationData;

    [[nodiscard]] unsigned int GetNumRows() const final {
        return column_data_[0].GetProbingTable().size();
    }
    [[nodiscard]] std::vector<int> GetTuple(int tuple_index) const;

    static std::unique_ptr<ColumnLayoutRelationData> CreateFrom(CSVParser& file_input,
                                                                bool is_null_eq_null,
                                                                int max_cols = -1,
                                                                long max_rows = -1);
};

