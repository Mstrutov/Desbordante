//
// Created by alexandrsmirn
//

#include "PruningMap.h"

PruningMap::PruningMap(RelationalSchema const* schema) {
    for (auto const& column : schema->getColumns()) {
        this->insert(std::make_pair(Vertical(*column), std::unordered_set<Vertical>()));
    }
}

void PruningMap::rebalance() {
    bool rebalancedGroup = false;

    do {
        rebalancedGroup = false;
        for (auto const& mapRow : *this) {
            Vertical const& key = mapRow.first;
            if (mapRow.second.size() > 1000) {
                rebalanceGroup(key);
                rebalancedGroup = true;
            }
        }
    } while (rebalancedGroup);
}

void PruningMap::rebalanceGroup(Vertical const& key) {
    std::unordered_set<Vertical> const& depsOfGroup = this->at(key);
    boost::dynamic_bitset<> invertedColumns = key.getColumnIndices().operator~();
    for (size_t columnIndex = invertedColumns.find_first();
         columnIndex < invertedColumns.size();
         columnIndex = invertedColumns.find_next(columnIndex))
    {
        //может быть не очень быстро?
        Vertical newKey = key.Union(*key.getSchema()->getColumn(columnIndex));
        std::unordered_set<Vertical> newGroup;
        this->insert(std::make_pair(newKey, newGroup));

        for (Vertical const& depOfGroup : depsOfGroup) {
            if (depOfGroup.contains(newKey)) {
                newGroup.insert(depOfGroup);
            }
        }
    }
    this->erase(key);
}

