#include "PliBasedFDAlgorithm.h"

void PliBasedFDAlgorithm::initialize() {
    relation_ = ColumnLayoutRelationData::createFrom(inputGenerator_, is_null_equal_null_);

    if (relation_->getColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");
    }
}
