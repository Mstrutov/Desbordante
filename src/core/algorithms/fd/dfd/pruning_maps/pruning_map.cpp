#include "pruning_map.h"

PruningMap::PruningMap(RelationalSchema const* schema) {
    for (auto const& column : schema->GetColumns()) {
        this->insert(std::make_pair(Vertical(column), std::unordered_set<Vertical>()));
    }
}

void PruningMap::Rebalance() {
    bool rebalanced_group = false;

    do {
        rebalanced_group = false;
        for (auto iter = this->begin(); iter != this->end();) {
            Vertical const& key = iter->first;
            auto const& related_verticals = iter->second;

            // RebalanceGroup() invalidates this iterator, because it erases the key element
            ++iter;
            if (related_verticals.size() > 1000) {
                RebalanceGroup(key);
                rebalanced_group = true;
            }
        }
    } while (rebalanced_group);
}

void PruningMap::RebalanceGroup(Vertical const& key) {
    auto const& deps_of_group = this->at(key);
    auto inverted_columns = key.GetColumnIndices().operator~();

    for (size_t column_index = inverted_columns.find_first();
         column_index < inverted_columns.size();
         column_index = inverted_columns.find_next(column_index)) {
        Vertical new_key = key.Union(*key.GetSchema()->GetColumn(column_index));
        std::unordered_set<Vertical> new_group;

        for (auto const& dep_of_group : deps_of_group) {
            if (dep_of_group.Contains(new_key)) {
                new_group.insert(dep_of_group);
            }
        }

        this->insert(std::make_pair(std::move(new_key), std::move(new_group)));
    }
    this->erase(key);
}
