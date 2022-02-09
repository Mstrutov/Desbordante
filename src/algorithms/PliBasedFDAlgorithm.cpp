#include "PliBasedFDAlgorithm.h"

void PliBasedFDAlgorithm::Initialize() {
    relation_ = ColumnLayoutRelationData::CreateFrom(input_generator_, is_null_equal_null_);

    if (relation_->GetColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");
    }
}

std::vector<Column const*> PliBasedFDAlgorithm::GetKeys() const {
    assert(relation_ != nullptr);

    std::vector<Column const*> keys;
    for (ColumnData const& col : relation_->GetColumnData()) {
        if (col.GetPositionListIndex()->GetNumNonSingletonCluster() == 0) {
            keys.push_back(col.GetColumn());
        }
    }

    return keys;
}
