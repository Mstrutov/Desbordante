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
    //std::variant<std::unique_ptr<PositionListIndex>, PositionListIndex*> positionListIndex_;
    std::shared_ptr<PositionListIndex> positionListIndex_;

public:
    ColumnData(Column const* column, std::unique_ptr<PositionListIndex> positionListIndex);
    std::vector<int> const& getProbingTable() const;
    int getProbingTableValue(int tupleIndex) const;
    Column const* getColumn() const { return column; }
    PositionListIndex const* getPositionListIndex() const;

    // Transfers positionListIndex_ ownership to the outside world. BE CAREFUL - other methods
    // of ColumnData get invalidated while the PLI is moved out
    // std::unique_ptr<PositionListIndex> moveOutPositionListIndex();

    std::shared_ptr<PositionListIndex> getPLIOwnership() { return positionListIndex_; }

    // Moves a PLI under the ownership of ColumnData
    // void moveInPositionListIndex(std::unique_ptr<PositionListIndex> positionListIndex ) { positionListIndex_ = std::move(positionListIndex); }

    //void shuffle();

    std::string toString() { return "Data for " + column->toString(); }
    bool operator==(const ColumnData& rhs);
};