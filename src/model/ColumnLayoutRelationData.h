//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <cmath>
#include <vector>

#include "ColumnData.h"
#include "RelationalSchema.h"
#include "RelationData.h"
#include "CSVParser.h"


using std::vector;

class ColumnLayoutRelationData : public RelationData {
private:
    vector<shared_ptr<ColumnData>> columnData;

public:
    vector<shared_ptr<ColumnData>> getColumnData() const override;
    shared_ptr<ColumnData> getColumnData(int columnIndex) const override;
    unsigned int getNumRows() const override;
    vector<int> getTuple(int tupleIndex) const override;
    void shuffleColumns() override;
    double getMaximumEntropy() { return std::log(getNumRows()); }

    ColumnLayoutRelationData(shared_ptr<RelationalSchema> const& schema, vector<shared_ptr<ColumnData>> columnData);

    static shared_ptr<ColumnLayoutRelationData> createFrom(CSVParser& fileInput, bool isNullEqNull);
    static shared_ptr<ColumnLayoutRelationData> createFrom(CSVParser& fileInput, bool isNullEqNull, int maxCols, long maxRows);
};