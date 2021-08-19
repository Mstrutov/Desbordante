#pragma once

#include "PositionListIndex.h"

class OrderedPartition {
private:
    //PositionListIndex const* const partition;
    double distinctiveness;
    //unsigned int numberOfRows;
    unsigned int columnIndex;
public:
    OrderedPartition() = delete;
    OrderedPartition(OrderedPartition const& other) = delete;
    OrderedPartition(OrderedPartition && other) = delete;
    OrderedPartition& operator=(OrderedPartition const& other) = delete;
    OrderedPartition& operator=(OrderedPartition && other) = delete;

    OrderedPartition(PositionListIndex const* const, unsigned int numberOfRows, unsigned int columnIndex);

    bool operator<(OrderedPartition const& other) const;
    unsigned int getColumnIndex() const { return this->columnIndex; }
};
