//
// Created by alex on 20.04.2021.
//

#include "NonDependenciesMap.h"
//#include "CustomComparator.h"

NonDependenciesMap::NonDependenciesMap(RelationalSchema const* schema) {
    for (auto const& column : schema->getColumns()) {
        this->insert(std::make_pair(Vertical(*column), vertical_set()));
    }
}

std::unordered_set<Vertical> NonDependenciesMap::getPrunedSupersets(std::unordered_set<Vertical> const& supersets) const {
    vertical_set prunedSupersets;
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
}

void NonDependenciesMap::addNewNonDependency(Vertical const& nodeToAdd) {
    for (auto const& mapRow : *this) {
        Vertical const& key = mapRow.first;

        if (nodeToAdd.contains(key)) {
            vertical_set nonDepsForKey = mapRow.second;
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
    //rebalance();
}

/*vector<shared_ptr<Vertical>>
NonDependenciesMap::getUncheckedSupersets(shared_ptr<Vertical> node, LatticeObservations & observations) {
    std::vector<shared_ptr<Vertical>> uncheckedSupersets;
    //dynamic_bitset<> invertedColumnIndices = node->getColumnIndices().flip();

    for (auto& subsetNode : node->getParents()) {
        if (observations.find(*subsetNode) == observations.end()) {
            uncheckedSupersets.push_back(std::move(subsetNode));
        }
    }

    for (size_t index = invertedColumnIndices.find_first(); index < invertedColumnIndices.size(); index = invertedColumnIndices.find_next(index)) {
        invertedColumnIndices[index] = false; //убираем одну из колонок
        auto supersetNode = std::make_shared<Vertical>(node->getSchema(), invertedColumnIndices.flip());
        auto supersetVerticalIter = this->find(*supersetNode); //TODO !!!лучше переделать без второго flip а просто бежать циклом по нулям

        if (supersetVerticalIter == this->end()) {
            if (!canBePruned(*supersetNode)) {
                uncheckedSupersets.push_back(std::move(supersetNode));
            } else {
                observations[*supersetNode] = NodeCategory::dependency;
            }
        }

        invertedColumnIndices[index] = true; //возвращаем как было
    }

    return uncheckedSupersets;
}*/