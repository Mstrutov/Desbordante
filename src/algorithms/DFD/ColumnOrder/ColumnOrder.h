//
// Created by alex on 11.07.2021.
//

#pragma once

#include "Vertical.h"
#include "ColumnData.h"

class ColumnOrder {
private:
    std::vector<int> order;
public:
    explicit ColumnOrder(ColumnLayoutRelationData const* const relationData);

    ColumnOrder() = default;
    //ColumnOrder(ColumnOrder const& other) = delete;
    //ColumnOrder& operator=(ColumnOrder const& other) = delete;

    std::vector<int> getOrderHighDistinctCount(Vertical const& columns) const;
    std::vector<int> getOrderLowDistinctCount(Vertical const& columns) const;
};
