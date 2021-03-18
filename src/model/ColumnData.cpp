//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#include "ColumnData.h"

#include <algorithm>
#include <random>
#include <utility>


ColumnData::ColumnData(Column const* column, std::unique_ptr<PositionListIndex> positionListIndex):
        column(column),
        positionListIndex_(std::move(positionListIndex)) {
        positionListIndex->forceCacheProbingTable();
    }

/*void ColumnData::shuffle() {
    std::random_device rd;
    std::mt19937 random(rd());
    std::shuffle(probingTable.begin(), probingTable.end(), random);
}*/

bool ColumnData::operator==(const ColumnData &rhs) {
    if (this == &rhs) return true;
    return this->column == rhs.column;
}