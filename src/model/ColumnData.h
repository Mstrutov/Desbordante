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
    Column const* column_;
    //std::variant<std::unique_ptr<PositionListIndex>, PositionListIndex*> position_list_index_;
    std::shared_ptr<util::PositionListIndex> position_list_index_;

public:
    ColumnData(Column const* column, std::unique_ptr<util::PositionListIndex> position_list_index);
    // Инвариант: конструктором гарантируется, что в ColumnData.PLI есть закешированная ProbingTable
    std::vector<int> const& GetProbingTable() const {
        return *position_list_index_->GetCachedProbingTable();
    }
    Column const* GetColumn() const { return column_; }
    int GetProbingTableValue(int tuple_index) const {
        return (*position_list_index_->GetCachedProbingTable())[tuple_index];
    }
    util::PositionListIndex const* GetPositionListIndex() const {
        return position_list_index_.get();
    }
    // TODO: посмотреть, что будет с производительностью, если добавить указатель на PT прямо сюда
    // по идее, это должно оптимизироваться инлайнингом

    // Transfers position_list_index_ ownership to the outside world. BE CAREFUL - other methods
    // of ColumnData get invalidated while the PLI is moved out
    // std::unique_ptr<PositionListIndex> moveOutPositionListIndex();

    std::shared_ptr<util::PositionListIndex> GetPliOwnership() { return position_list_index_; }
    std::shared_ptr<util::PositionListIndex const> GetPliOwnership() const {
        return position_list_index_;
    }

    // Moves a PLI under the ownership of ColumnData
    // void moveInPositionListIndex(std::unique_ptr<PositionListIndex> positionListIndex ) { position_list_index_ = std::move(positionListIndex); }

    //void shuffle();

    std::string ToString() { return "Data for " + column_->ToString(); }
    bool operator==(const ColumnData& rhs);
};
