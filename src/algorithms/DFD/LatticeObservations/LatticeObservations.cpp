//
// Created by alex on 13.04.2021.
//

#include "LatticeObservations.h"

NodeCategory LatticeObservations::updateDependencyCategory(Vertical const& vertical) {
    //бежим по множеству подмножеств и смотрим. если получили все независимости, то это мин зависимость

    if (vertical.getArity() > 1) {
        boost::dynamic_bitset<> columnIndices = vertical.getColumnIndices(); //копируем индексы
        bool hasUncheckedSubset = false;

        for (size_t index = columnIndices.find_first();
             index < columnIndices.size(); index = columnIndices.find_next(index)) {
            columnIndices[index] = false; //убираем одну из колонок
            auto const subsetVerticalIter = this->find(
                    Vertical(vertical.getSchema(), columnIndices)); //TODO передаем временнй объект??

            //если какое-то подмножество не посещено либо оно не является антизависимостью, то не подходит

            if (subsetVerticalIter == this->end()) {
                //если нашли нерассмотренное подмножество
                hasUncheckedSubset = true;
            } else {
                NodeCategory const &subsetVerticalCategory = subsetVerticalIter->second;
                if (subsetVerticalCategory == NodeCategory::minimalDependency ||
                    subsetVerticalCategory == NodeCategory::dependency ||
                    subsetVerticalCategory == NodeCategory::candidateMinimalDependency
                ) {
                    return NodeCategory::dependency;
                }
            }
            //если все норм
            columnIndices[index] = true; //возвращаем как было
        }
        return hasUncheckedSubset ? NodeCategory::candidateMinimalDependency : NodeCategory::minimalDependency;
    } else {
        return NodeCategory::minimalDependency;
    }
}

NodeCategory LatticeObservations::updateNonDependencyCategory(Vertical const& vertical, int rhsIndex) {
    boost::dynamic_bitset<> columnIndices = vertical.getColumnIndices();
    //columnIndices[rhsIndex] = true;
    bool hasUncheckedSuperset = false;

    for (size_t index = 0; index < columnIndices.size(); index++) {
        if (!columnIndices[index] && index != rhsIndex) {
            columnIndices[index] = true; //убираем одну из колонок
            auto const supersetVerticalIter = this->find(Vertical(vertical.getSchema(), columnIndices)); //TODO !!!лучше переделать без второго flip а просто бежать циклом по нулям

            if (supersetVerticalIter == this->end()) {
                //если нашли нерассмотренное надмножество
                hasUncheckedSuperset = true;
            } else {
                NodeCategory const &supersetVerticalCategory = supersetVerticalIter->second;
                if (supersetVerticalCategory == NodeCategory::maximalNonDependency ||
                    supersetVerticalCategory == NodeCategory::nonDependency ||
                    supersetVerticalCategory == NodeCategory::candidateMaximalNonDependency
                ) {
                    return NodeCategory::nonDependency;
                }
            }
            columnIndices[index] = false; //возвращаем как было
        }
    }
    return hasUncheckedSuperset ? NodeCategory::candidateMaximalNonDependency : NodeCategory::maximalNonDependency;
}

bool LatticeObservations::isCandidate(Vertical const& node) {
    auto nodeIter = this->find(node);
    if (nodeIter == this->end()) {
        return false;
    } else {
        return nodeIter->second == NodeCategory::candidateMaximalNonDependency ||
               nodeIter->second == NodeCategory::candidateMinimalDependency;
    }
}


std::unordered_set<Vertical>
LatticeObservations::getUncheckedSubsets(Vertical const& node, size_t rhsIndex) {
    /*vertical_set uncheckedSubsets;

    for (auto& subsetNode : node->getParents()) {
        if (this->find(*subsetNode) == this->end()) {
            uncheckedSubsets.insert(std::move(subsetNode));
        }
    }*/

    boost::dynamic_bitset<> indices = node.getColumnIndices();
    std::unordered_set<Vertical> uncheckedSubsets;

    for (size_t index = 0; index < indices.size(); index++) {
        if (index != rhsIndex && indices[index]) {
            indices[index] = false;
            Vertical subsetNode = Vertical(node.getSchema(), indices);

            if (this->find(subsetNode) == this->end()) {
                uncheckedSubsets.insert(std::move(subsetNode));
            }
            indices[index] = true;
        }
    }

    return uncheckedSubsets;
}

std::unordered_set<Vertical>
LatticeObservations::getUncheckedSupersets(Vertical const& node, size_t rhsIndex) {
    boost::dynamic_bitset<> indices = node.getColumnIndices();
    std::unordered_set<Vertical> uncheckedSupersets;

    for (size_t index = 0; index < indices.size(); index++) {
        if (index != rhsIndex && !indices[index]) {
            indices[index] = true;
            Vertical supersetNode = Vertical(node.getSchema(), indices);

            if (this->find(supersetNode) == this->end()) {
                uncheckedSupersets.insert(std::move(supersetNode));
            }
            indices[index] = false;
        }
    }
    return uncheckedSupersets;
}


