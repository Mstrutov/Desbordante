//
// Created by alexandrsmirn
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

    std::vector<int> getOrderHighDistinctCount(Vertical const& columns) const;
    std::vector<int> getOrderLowDistinctCount(Vertical const& columns) const;
};
