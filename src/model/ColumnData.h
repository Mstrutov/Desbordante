//
// Created by Ilya Vologin
// https://github.com/cupertank
//

#pragma once

#include <vector>

#include "AbstractColumnData.h"
#include "Column.h"
#include "PositionListIndex.h"

class ColumnData final : public model::AbstractColumnData {
    std::shared_ptr<util::PositionListIndex> position_list_index_;

public:
    ColumnData(Column const* column, std::unique_ptr<util::PositionListIndex> position_list_index)
        : AbstractColumnData(column), position_list_index_(std::move(position_list_index)) {
        position_list_index_->ForceCacheProbingTable();
    }
    // Инвариант: конструктором гарантируется, что в ColumnData.PLI есть закешированная ProbingTable
    std::vector<int> const& GetProbingTable() const {
        return *position_list_index_->GetCachedProbingTable();
    }
    int GetProbingTableValue(int tuple_index) const {
        return (*position_list_index_->GetCachedProbingTable())[tuple_index];
    }
    util::PositionListIndex const* GetPositionListIndex() const {
        return position_list_index_.get();
    }

    std::shared_ptr<util::PositionListIndex> GetPliOwnership() { return position_list_index_; }
    std::shared_ptr<util::PositionListIndex const> GetPliOwnership() const {
        return position_list_index_;
    }

    std::string ToString() const final { return "Data for " + column_->ToString(); }
};

