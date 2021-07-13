//
// Created by alex on 11.07.2021.
//

#include "ColumnOrder.h"

#include "OrderedPartition.h"
#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"

ColumnOrder::ColumnOrder(ColumnLayoutRelationData const* const relationData)
        : order(relationData->getSchema()->getNumColumns()) { //TODO check
    std::set<OrderedPartition> partitions;
    for (auto const& columnData : relationData->getColumnData()) {
        partitions.emplace(columnData.getPositionListIndex(), relationData->getNumRows(), columnData.getColumn()->getIndex());
    }

    int orderIndex = 0;
    for (OrderedPartition const& partition : partitions) {
        order[orderIndex++] = partition.getColumnIndex();
    }
}

std::vector<int> ColumnOrder::getOrderHighDistinctCount(const Vertical &columns) const {
    std::vector<int> orderForColumns(columns.getArity()); //TODO возможно getArity прожорливый?

    int currentOrderIndex = 0;
    for (int i = 0; i < this->order.size(); ++i) {
        if (columns.getColumnIndices()[order[i]]) {
            orderForColumns[currentOrderIndex++] = this->order[i];
        }
    }

    return orderForColumns;
}

std::vector<int> ColumnOrder::getOrderLowDistinctCount(const Vertical &columns) const {
    std::vector<int> orderForColumns(columns.getArity()); //TODO возможно getArity прожорливый?

    int currentOrderIndex = 0;
    for (int i = this->order.size() - 1; i >= 0; --i) {
        if (columns.getColumnIndices()[order[i]]) {
            orderForColumns[currentOrderIndex++] = this->order[i];
        }
    }

    return orderForColumns;
}
