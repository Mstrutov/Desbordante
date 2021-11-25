//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#include "ColumnData.h"

#include <algorithm>
#include <random>
#include <utility>


ColumnData::ColumnData(Column const* column,
                       std::unique_ptr<util::PositionListIndex> positionListIndex):
        column(column),
        positionListIndex_(std::move(positionListIndex)) {
        //std::get<std::unique_ptr<PositionListIndex>>(positionListIndex_)->forceCacheProbingTable();
        positionListIndex_->forceCacheProbingTable();
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

//std::vector<int> const &ColumnData::getProbingTable() const {
//    /*if (std::holds_alternative<std::unique_ptr<PositionListIndex>>(positionListIndex_)) {
//        return *std::get<std::unique_ptr<PositionListIndex>>(positionListIndex_)->getProbingTable();
//    } else {
//        return *std::get<PositionListIndex*>(positionListIndex_)->getProbingTable();
//    }*/
//    return *positionListIndex_->getProbingTable();
//}
//
//int ColumnData::getProbingTableValue(int tupleIndex) const {
////    if (std::holds_alternative<std::unique_ptr<PositionListIndex>>(positionListIndex_)) {
////        return (*std::get<std::unique_ptr<PositionListIndex>>(positionListIndex_)->getProbingTable())[tupleIndex];
////    } else {
////        return (*std::get<PositionListIndex*>(positionListIndex_)->getProbingTable())[tupleIndex];
////    }
//
//}
//
//PositionListIndex const *ColumnData::getPositionListIndex() const {
////    if (std::holds_alternative<std::unique_ptr<PositionListIndex>>(positionListIndex_)) {
////        return std::get<std::unique_ptr<PositionListIndex>>(positionListIndex_).get();
////    } else {
////        return std::get<PositionListIndex*>(positionListIndex_);
////    }
//    return positionListIndex_.get();
//}

//std::unique_ptr<PositionListIndex> ColumnData::moveOutPositionListIndex() {
//    if (std::holds_alternative<std::unique_ptr<PositionListIndex>>(positionListIndex_)) {
//        auto ownedPtr = std::move(std::get<std::unique_ptr<PositionListIndex>>(positionListIndex_));
//        positionListIndex_ = ownedPtr.get();
//        return ownedPtr;
//    } else {
//        throw std::logic_error("PLI is already moved out");
//    }
//}
