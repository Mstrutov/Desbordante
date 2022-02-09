#include "OrderedPartition.h"

OrderedPartition::OrderedPartition(util::PositionListIndex const* const partition,
                                   unsigned int number_of_rows, unsigned int column_index)
    : column_index_(column_index) {
    this->distinctiveness_ =
        (double)(number_of_rows - partition->GetNumNonSingletonCluster()) / number_of_rows;
}

bool OrderedPartition::operator<(OrderedPartition const& other) const {
    if (this->distinctiveness_ == other.distinctiveness_) {
        return this->column_index_ < other.column_index_;
    }
    return this->distinctiveness_ < other.distinctiveness_;
}
