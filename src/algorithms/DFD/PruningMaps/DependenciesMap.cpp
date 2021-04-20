//
// Created by alexandrsmirn
//

#include "DependenciesMap.h"

DependenciesMap::DependenciesMap(shared_ptr<RelationalSchema> schema) {
    for (auto const& column : schema->getColumns()) {
        this->insert(std::make_pair(Vertical(*column), std::unordered_set<shared_ptr<Vertical>>()));
    }
}

void DependenciesMap::addNewDependency(shared_ptr<Vertical> node) { //пока версия без балансировок
    using std::unordered_set;

    for (auto const& column : node->getColumns()) {
        unordered_set<shared_ptr<Vertical>>& verticalSet = this->find(Vertical(*column))->second;
        for (auto const& vertical : verticalSet) {
            if (vertical->contains(*node)) {
                verticalSet.erase(vertical);    //удаляем надмножества
            }
        }
        verticalSet.insert(node);
    }
}

std::vector<shared_ptr<Vertical>> DependenciesMap::getUncheckedSubsets(shared_ptr<Vertical> node, LatticeObservations const& observations) {
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

bool DependenciesMap::canBePruned(Vertical const& node) {
    using std::unordered_set;

    for (auto const& column : node.getColumns()) {
        unordered_set<shared_ptr<Vertical>>& verticalSet = this->find(Vertical(*column))->second;
        for (auto const& vertical : verticalSet) {
            if (node.contains(*vertical)) { //TODO contains проверяет на равенство?
                return true;
            }
        }
    }

    return false;
}


