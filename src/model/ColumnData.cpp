//
// Created by kek on 12.07.19.
//

#include "ColumnData.h"

#include <algorithm>
#include <random>
#include <utility>

using namespace std;

ColumnData::ColumnData(shared_ptr<Column>& column, vector<int> probingTable, shared_ptr<PositionListIndex>& positionListIndex):
    column(column),
    probingTable(std::make_unique<vector<int>>(std::move(probingTable))),
    positionListIndex(positionListIndex)
    {}

//TODO: Random проверь
void ColumnData::shuffle() {
    random_device rd;
    mt19937 random(rd());
    std::shuffle(probingTable->begin(), probingTable->end(), random);
}


bool ColumnData::operator==(const ColumnData &rhs) {
    if (this == &rhs) return true;
    return this->column == rhs.column;
}