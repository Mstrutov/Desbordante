//
// Created by Ilya Vologin
// https://github.com/cupertank
//


#include "ColumnData.h"

#include <algorithm>
#include <random>
#include <utility>

ColumnData::ColumnData(Column const* column,
                       std::unique_ptr<util::PositionListIndex> position_list_index)
    : column_(column), position_list_index_(std::move(position_list_index)) {
    // std::get<std::unique_ptr<PositionListIndex>>(position_list_index_)->ForceCacheProbingTable();
    position_list_index_->ForceCacheProbingTable();
}

/*void ColumnData::shuffle() {
    std::random_device rd;
    std::mt19937 random(rd());
    std::shuffle(probingTable.begin(), probingTable.end(), random);
}*/

bool ColumnData::operator==(const ColumnData& rhs) {
    if (this == &rhs) return true;
    return this->column_ == rhs.column_;
}

//std::vector<int> const &ColumnData::GetProbingTable() const {
//    /*if (std::holds_alternative<std::unique_ptr<PositionListIndex>>(position_list_index_)) {
//        return *std::get<std::unique_ptr<PositionListIndex>>(position_list_index_)->GetProbingTable();
//    } else {
//        return *std::get<PositionListIndex*>(position_list_index_)->GetProbingTable();
//    }*/
//    return *position_list_index_->GetProbingTable();
//}
//
//int ColumnData::getProbingTableValue(int tupleIndex) const {
////    if (std::holds_alternative<std::unique_ptr<PositionListIndex>>(position_list_index_)) {
////        return (*std::get<std::unique_ptr<PositionListIndex>>(position_list_index_)->GetProbingTable())[tupleIndex];
////    } else {
////        return (*std::get<PositionListIndex*>(position_list_index_)->GetProbingTable())[tupleIndex];
////    }
//
//}
//
//PositionListIndex const *ColumnData::GetPositionListIndex() const {
////    if (std::holds_alternative<std::unique_ptr<PositionListIndex>>(position_list_index_)) {
////        return std::get<std::unique_ptr<PositionListIndex>>(position_list_index_).get();
////    } else {
////        return std::get<PositionListIndex*>(position_list_index_);
////    }
//    return position_list_index_.get();
//}

//std::unique_ptr<PositionListIndex> ColumnData::moveOutPositionListIndex() {
//    if (std::holds_alternative<std::unique_ptr<PositionListIndex>>(position_list_index_)) {
//        auto ownedPtr = std::move(std::get<std::unique_ptr<PositionListIndex>>(position_list_index_));
//        position_list_index_ = ownedPtr.get();
//        return ownedPtr;
//    } else {
//        throw std::logic_error("PLI is already moved out");
//    }
//}
