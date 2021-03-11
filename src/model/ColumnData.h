//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <vector>

#include "Column.h"
#include "PositionListIndex.h"

class ColumnData {
private:
    Column const* column;
    std::unique_ptr<PositionListIndex> positionListIndex;

public:
    ColumnData(Column const* column, std::unique_ptr<PositionListIndex> positionListIndex);
    std::vector<int> const& getProbingTable() const { return (*positionListIndex->getProbingTable()); }
    int getProbingTableValue(int tupleIndex) const { return (*positionListIndex->getProbingTable())[tupleIndex]; }
    Column const* getColumn() const { return column; }
    PositionListIndex const* getPositionListIndex() const { return positionListIndex.get(); }

    // Transfers positionListIndex ownership to the outside world. BE CAREFUL - other methods
    // of ColumnData get invalidated while the PLI is moved out
    // std::unique_ptr<PositionListIndex> movePositionListIndex() { return std::move(positionListIndex); }

    //void shuffle();

    std::string toString() { return "Data for " + column->toString(); }
    bool operator==(const ColumnData& rhs);
};