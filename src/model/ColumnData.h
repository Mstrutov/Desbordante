//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#pragma once

#include <vector>

#include "Column.h"
#include "PositionListIndex.h"

using std::vector;

class ColumnData {
private:
    shared_ptr<Column> column;
    std::unique_ptr<vector<int>> probingTable;
    shared_ptr<PositionListIndex> positionListIndex;

public:
    ColumnData(shared_ptr<Column>& column, vector<int> probingTable, shared_ptr<PositionListIndex>& positionListIndex);
    vector<int>* getProbingTable() { return probingTable.get(); }
    int getProbingTableValue(int tupleIndex) { return (*probingTable)[tupleIndex]; }
    shared_ptr<Column> getColumn() { return column; }
    shared_ptr<PositionListIndex> getPositionListIndex() { return positionListIndex; }
    void shuffle();
    string toString() { return "Data for " + column->toString(); }
    bool operator==(const ColumnData& rhs);
};