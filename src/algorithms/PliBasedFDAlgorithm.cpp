#include "PliBasedFDAlgorithm.h"

void PliBasedFDAlgorithm::initialize() {
    relation_ = ColumnLayoutRelationData::createFrom(inputGenerator_, is_null_equal_null_);

    if (relation_->getColumnData().empty()) {
        throw std::runtime_error("Got an empty .csv file: FD mining is meaningless.");
    }
}

std::vector<Column const*> PliBasedFDAlgorithm::getKeys() const {
   assert(relation_ != nullptr);

   std::vector<Column const*> keys;
   for (ColumnData const& col : relation_->getColumnData()) {
       if (col.getPositionListIndex()->getNumNonSingletonCluster() == 0) {
           keys.push_back(col.getColumn());
       }
   }

   return keys;
}
