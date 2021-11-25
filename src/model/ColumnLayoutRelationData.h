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
    std::vector<ColumnData> columnData;

public:
    std::vector<ColumnData>& getColumnData() override { return columnData; };
    std::vector<ColumnData> const& getColumnData() const override { return columnData; };
    ColumnData& getColumnData(int columnIndex) override;
    ColumnData const& getColumnData(int columnIndex) const override {
        return columnData[columnIndex];
    }
    unsigned int getNumRows() const override { return columnData[0].getProbingTable().size(); }
    std::vector<int> getTuple(int tupleIndex) const override;

    //void shuffleColumns() override = 0;

    double getMaximumEntropy() const { return std::log(getNumRows()); }

    ColumnLayoutRelationData(std::unique_ptr<RelationalSchema> schema,
                             std::vector<ColumnData> columnData);

    static std::unique_ptr<ColumnLayoutRelationData> createFrom(CSVParser& fileInput,
                                                                bool isNullEqNull);
    static std::unique_ptr<ColumnLayoutRelationData> createFrom(
            CSVParser& fileInput, bool isNullEqNull, int maxCols, long maxRows);
};
