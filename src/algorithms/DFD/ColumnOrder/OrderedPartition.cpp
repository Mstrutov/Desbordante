#include "OrderedPartition.h"

OrderedPartition::OrderedPartition(const PositionListIndex *const partition, unsigned int numberOfRows, unsigned int columnIndex)
        : PositionListIndex(*partition), numberOfRows(numberOfRows), columnIndex(columnIndex) {
    this->distinctiveness = (double)(this->numberOfRows - this->getNumNonSingletonCluster()) / this->numberOfRows;
}

bool OrderedPartition::operator<(OrderedPartition const& other) const {
    if (this->distinctiveness == other.distinctiveness) {
        return this->columnIndex < other.columnIndex;
    }
    return this->distinctiveness < other.distinctiveness;
}
