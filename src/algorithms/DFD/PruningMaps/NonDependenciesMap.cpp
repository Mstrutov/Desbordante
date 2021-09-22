#include "NonDependenciesMap.h"

NonDependenciesMap::NonDependenciesMap(RelationalSchema const* schema)
    : PruningMap(schema) {}

std::unordered_set<Vertical> NonDependenciesMap::getPrunedSupersets(std::unordered_set<Vertical> const& supersets) const {
    std::unordered_set<Vertical> prunedSupersets;
    for (auto const& node : supersets) {
        if (canBePruned(node)) {
            prunedSupersets.insert(node);
        }
    }
    return prunedSupersets;
}

bool NonDependenciesMap::canBePruned(const Vertical &node) const {
    for (auto const& mapRow : *this) {
        Vertical const& key = mapRow.first;
        if (node.contains(key)) {
            for (Vertical const& nonDependency : mapRow.second) {
                if (nonDependency.contains(node)) {
                    return true;
                }
            }
        }
    }
    return false;
}

void NonDependenciesMap::addNewNonDependency(Vertical const& nodeToAdd) {
    for (auto& mapRow : *this) {
        Vertical const& key = mapRow.first;

        if (nodeToAdd.contains(key)) {
            auto& nonDepsForKey = mapRow.second;
            bool hasSupersetEntry = false;

            for (auto iter = nonDepsForKey.begin(); iter != nonDepsForKey.end(); ) {
                //если совпадают, то contains = true
                Vertical const& nonDep = *iter;
                if (nonDep.contains(nodeToAdd)) {
                    hasSupersetEntry = true;
                    break;
                } else if (nodeToAdd.contains(nonDep)) {
                    iter = nonDepsForKey.erase(iter);
                } else {
                    iter++;
                }
            }

            if (!hasSupersetEntry) {
                nonDepsForKey.insert(nodeToAdd);
            }
        }
    }
    rebalance();
}
