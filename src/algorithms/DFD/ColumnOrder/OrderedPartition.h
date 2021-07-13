//
// Created by alex on 11.07.2021.
//

#pragma once

#include "PositionListIndex.h"

class OrderedPartition : public PositionListIndex {
private:
    double distinctiveness;
    unsigned int numberOfRows;
    unsigned int columnIndex;
public:
    OrderedPartition() = delete;
    OrderedPartition(OrderedPartition const& other) = delete;
    OrderedPartition& operator=(OrderedPartition const& other) = delete;

    OrderedPartition(PositionListIndex const* const, unsigned int numberOfRows, unsigned int columnIndex);
    bool operator<(OrderedPartition const& other) const;
    unsigned int getColumnIndex() const { return this->columnIndex; }
};
