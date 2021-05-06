//
// Created by alexandrsmirn
//

#include "DependenciesMap.h"

DependenciesMap::DependenciesMap(shared_ptr<RelationalSchema> schema) {
    for (auto const& column : schema->getColumns()) {
        this->insert(std::make_pair(Vertical(*column), vertical_set()));
    }
}

DependenciesMap::vertical_set DependenciesMap::getPrunedSubsets(vertical_set const& subsets) const {
    vertical_set prunedSubsets;
    for (auto const& node : subsets) {
        if (canBePruned(*node)) {
            prunedSubsets.insert(node);
        }
    }
    return prunedSubsets;
}

void DependenciesMap::addNewDependency(shared_ptr<Vertical> nodeToAdd) { //пока версия без балансировок
    for (auto const& mapRow : *this) {
        Vertical const& key = mapRow.first;

        if (nodeToAdd->contains(key)) {
            vertical_set depsForKey = mapRow.second;
            bool hasSubsetEntry = false;

            for (auto iter = depsForKey.begin(); iter != depsForKey.end(); ) {
                //если совпадают, то contains = true
                shared_ptr<Vertical> const& dep = *iter;
                if (nodeToAdd->contains(*dep)) {
                    hasSubsetEntry = true;
                    break;
                } else if (dep->contains(*nodeToAdd)) {
                    iter = depsForKey.erase(iter);
                } else {
                    iter++;
                }
            }

            if (!hasSubsetEntry) {
                depsForKey.insert(nodeToAdd);
            }
        }
    }
    //rebalance();
}

/*std::vector<shared_ptr<Vertical>> DependenciesMap::getUncheckedSubsets(shared_ptr<Vertical> node, LatticeObservations & observations) {
    std::vector<shared_ptr<Vertical>> uncheckedSubsets;

    for (auto& subsetNode : node->getParents()) {
        if (observations.find(*subsetNode) == observations.end()) {
            if (!canBePruned(*subsetNode)) {
                uncheckedSubsets.push_back(std::move(subsetNode));
            } else {
                observations[*subsetNode] = NodeCategory::nonDependency;
            }
        }
    }

    return uncheckedSubsets;
}*/

bool DependenciesMap::canBePruned(Vertical const& node) const {
    for (auto const& mapRow : *this) {
        Vertical const& key = mapRow.first;
        if (node.contains(key)) {
            for (shared_ptr<Vertical> const& dependency : mapRow.second) {
                if (node.contains(*dependency)) {
                    return true;
                }
            }
        }
    }
}
