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
        for (auto iter = this->begin(); iter != this->end(); ) {
            Vertical const& key = iter->first;
            auto const& relatedVerticals = iter->second;

            //метод rebalanceGroup инвалидирует этот итератор, т.к. удаляет элемент
            ++iter;
            if (relatedVerticals.size() > 1000) {
                rebalanceGroup(key);
                rebalancedGroup = true;
            }
        }
    } while (rebalancedGroup);
}

void PruningMap::rebalanceGroup(Vertical const& key) {
    auto const& depsOfGroup = this->at(key);
    auto invertedColumns = key.getColumnIndices().operator~();

    for (size_t columnIndex = invertedColumns.find_first();
         columnIndex < invertedColumns.size();
         columnIndex = invertedColumns.find_next(columnIndex))
    {
        Vertical newKey = key.Union(*key.getSchema()->getColumn(columnIndex));
        std::unordered_set<Vertical> newGroup;

        for (auto const& depOfGroup : depsOfGroup) {
            if (depOfGroup.contains(newKey)) {
                newGroup.insert(depOfGroup);
            }
        }

        this->insert(std::make_pair(std::move(newKey), std::move(newGroup)));
    }
    this->erase(key);
}
