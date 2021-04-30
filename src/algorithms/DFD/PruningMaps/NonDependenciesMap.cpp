//
// Created by alex on 20.04.2021.
//

#include "NonDependenciesMap.h"

NonDependenciesMap::NonDependenciesMap(shared_ptr<RelationalSchema> schema) {
    for (auto const& column : schema->getColumns()) {
        this->insert(std::make_pair(Vertical(*column), std::unordered_set<shared_ptr<Vertical>>()));
    }
}

void NonDependenciesMap::addNewNonDependency(shared_ptr<Vertical> node) {
    using std::unordered_set;

    for (auto const& column : node->getColumns()) {
        unordered_set<shared_ptr<Vertical>>& verticalSet = this->find(Vertical(*column))->second;
        /*for (auto const& vertical : verticalSet) {
            if (node->contains(*vertical)) {
                verticalSet.erase(vertical);    //удаляем подмножества
            }
        }*/

        for (auto iter = verticalSet.begin(); iter != verticalSet.end(); ) {
            if (node->contains(**iter)) {
                iter = verticalSet.erase(iter); //удаляем подмножества
            } else {
                iter++;
            }
        }

        verticalSet.insert(node);
    }
}

vector<shared_ptr<Vertical>>
NonDependenciesMap::getUncheckedSupersets(shared_ptr<Vertical> node, LatticeObservations const& observations) {
    std::vector<shared_ptr<Vertical>> uncheckedSubsets;

    for (auto& subsetNode : node->getParents()) {
        if (observations.find(*subsetNode) == observations.end() &&
            !canBePruned(*subsetNode)
        ) {
            uncheckedSubsets.push_back(std::move(subsetNode));
        }
    }

    return uncheckedSubsets;
}

bool NonDependenciesMap::canBePruned(const Vertical &node) {
    using std::unordered_set;

    for (auto const& column : node.getColumns()) {
        unordered_set<shared_ptr<Vertical>>& verticalContainingColumnSet = this->find(Vertical(*column))->second;
        for (auto const& vertical : verticalContainingColumnSet) {
            if (vertical->contains(node)) { //TODO contains проверяет на равенство?
                return true;
            }
        }
    }
    return false;
}
