#include "ColumnOrder.h"

#include "OrderedPartition.h"
#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"

ColumnOrder::ColumnOrder(ColumnLayoutRelationData const* const relation_data)
    : order_(relation_data->GetSchema()->GetNumColumns()) {
    std::set<OrderedPartition> partitions;
    for (auto const& column_data : relation_data->GetColumnData()) {
        partitions.emplace(column_data.GetPositionListIndex(), relation_data->GetNumRows(),
                           column_data.GetColumn()->GetIndex());
    }

    int order_index = 0;
    for (auto const& partition : partitions) {
        order_[order_index++] = partition.GetColumnIndex();
    }
}

std::vector<int> ColumnOrder::GetOrderHighDistinctCount(const Vertical& columns) const {
    std::vector<int> order_for_columns(columns.GetArity());

    int current_order_index = 0;
    for (int column_index : order_) {
        if (columns.GetColumnIndices()[column_index]) {
            order_for_columns[current_order_index++] = column_index;
        }
    }

    return order_for_columns;
}

std::vector<int> ColumnOrder::GetOrderLowDistinctCount(const Vertical& columns) const {
    std::vector<int> order_for_columns(columns.GetArity());

    assert(!order_.empty());
    int current_order_index = 0;
    for (int i = this->order_.size() - 1; i >= 0; --i) {
        if (columns.GetColumnIndices()[order_[i]]) {
            order_for_columns[current_order_index++] = this->order_[i];
        }
    }

    return order_for_columns;
}
