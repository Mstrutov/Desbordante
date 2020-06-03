//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <cmath>
#include <vector>

#include "model/ColumnData.h"
#include "model/RelationalSchema.h"
#include "model/RelationData.h"
#include "parser/CSVParser.h"


using std::vector;

class ColumnLayoutRelationData : public RelationData {
private:
    vector<shared_ptr<ColumnData>> columnData;
    ColumnLayoutRelationData(shared_ptr<RelationalSchema>& schema, vector<shared_ptr<ColumnData>> columnData);

public:
    vector<shared_ptr<ColumnData>> getColumnData() override;
    shared_ptr<ColumnData> getColumnData(int columnIndex) override;
    unsigned int getNumRows() override;
    vector<int> getTuple(int tupleIndex) override;
    void shuffleColumns() override;
    double getMaximumEntropy() { return std::log(getNumRows()); }

    static shared_ptr<ColumnLayoutRelationData> createFrom(CSVParser& fileInput, bool isNullEqNull);
    static shared_ptr<ColumnLayoutRelationData> createFrom(CSVParser& fileInput, bool isNullEqNull, int maxCols, long maxRows);
};