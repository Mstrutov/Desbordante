#include "OrderedPartition.h"

OrderedPartition::OrderedPartition(util::PositionListIndex const *const partition,
                                   unsigned int numberOfRows, unsigned int columnIndex)
        : columnIndex(columnIndex) {
    this->distinctiveness = (double)(numberOfRows - partition->getNumNonSingletonCluster()) / numberOfRows;
}

bool OrderedPartition::operator<(OrderedPartition const& other) const {
    if (this->distinctiveness == other.distinctiveness) {
        return this->columnIndex < other.columnIndex;
    }
    return this->distinctiveness < other.distinctiveness;
}
