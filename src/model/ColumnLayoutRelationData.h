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


class ColumnLayoutRelationData : public RelationData {
private:
    std::vector<ColumnData> column_data_;

public:
    std::vector<ColumnData>& GetColumnData() override { return column_data_; };
    std::vector<ColumnData> const& GetColumnData() const override { return column_data_; };
    ColumnData& GetColumnData(int column_index) override;
    ColumnData const& GetColumnData(int column_index) const override {
        return column_data_[column_index];
    }
    unsigned int GetNumRows() const override { return column_data_[0].GetProbingTable().size(); }
    std::vector<int> GetTuple(int tuple_index) const override;

    //void shuffleColumns() override = 0;

    double GetMaximumEntropy() const { return std::log(GetNumRows()); }

    ColumnLayoutRelationData(std::unique_ptr<RelationalSchema> schema,
                             std::vector<ColumnData> column_data);

    static std::unique_ptr<ColumnLayoutRelationData> CreateFrom(CSVParser& file_input,
                                                                bool is_null_eq_null);
    static std::unique_ptr<ColumnLayoutRelationData> CreateFrom(
        CSVParser& file_input, bool is_null_eq_null, int max_cols, long max_rows);
};
