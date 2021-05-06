//
// Created by alex on 13.04.2021.
//

#include "LatticeObservations.h"

NodeCategory LatticeObservations::updateDependencyCategory(const shared_ptr<Vertical> &vertical) {
    //бежим по множеству подмножеств и смотрим. если получили все независимости, то это мин зависимость
    //TODO что если узел единичный
    dynamic_bitset<> columnIndices = vertical->getColumnIndices(); //копируем индексы
    bool hasUncheckedSubset = false;

    for (size_t index = columnIndices.find_first(); index < columnIndices.size(); index = columnIndices.find_next(index)) {
        columnIndices[index] = false; //убираем одну из колонок
        auto const subsetVerticalIter = this->find(Vertical(vertical->getSchema(), columnIndices)); //TODO передаем временнй объект??

        //если какое-то подмножество не посещено либо оно не является антизависимостью, то не подходит
        //TODO может быть такое говно что если нашли какую-то неантизависимость, то убираем кандидатство?

        if (subsetVerticalIter == this->end()) {
            //если нашли нерассмотренное подмножество
            hasUncheckedSubset = true;
        } else {
            NodeCategory const& subsetVerticalCategory = subsetVerticalIter->second;
            if (subsetVerticalCategory == NodeCategory::minimalDependency ||
                subsetVerticalCategory == NodeCategory::dependency
            ) {
                return NodeCategory::dependency;
            }
        }
        //если все норм
        columnIndices[index] = true; //возвращаем как было
    }
    return hasUncheckedSubset ? NodeCategory::candidateMaximalNonDependency : NodeCategory::minimalDependency;
}

NodeCategory LatticeObservations::updateNonDependencyCategory(const shared_ptr<Vertical> &vertical) {
    //копируем индексы колонок, которые не содержатся в искомом вертикале, чтобы исследовать надмножетва
    dynamic_bitset<> invertedColumnIndices = vertical->getColumnIndices().flip();
    bool hasUncheckedSuperset = false;

    for (size_t index = invertedColumnIndices.find_first(); index < invertedColumnIndices.size(); index = invertedColumnIndices.find_next(index)) {
        invertedColumnIndices[index] = false; //убираем одну из колонок
        auto const supersetVerticalIter = this->find(Vertical(vertical->getSchema(), invertedColumnIndices.flip())); //TODO !!!лучше переделать без второго flip а просто бежать циклом по нулям

        if (supersetVerticalIter == this->end()) {
            //если нашли нерассмотренное надмножество
            hasUncheckedSuperset = true;
        } else {
            NodeCategory const& supersetVerticalCategory = supersetVerticalIter->second;
            if (supersetVerticalCategory == NodeCategory::maximalNonDependency ||
                supersetVerticalCategory == NodeCategory::nonDependency
            ) {
                return NodeCategory::nonDependency;
            }
        }
        //если все норм
        invertedColumnIndices[index] = true; //возвращаем как было
    }
    return hasUncheckedSuperset ? NodeCategory::candidateMaximalNonDependency : NodeCategory::maximalNonDependency;
}

bool LatticeObservations::inferCategory(shared_ptr<Vertical> const& node) { //TODO можно ли оптимизировать?
    /*dynamic_bitset<> columnIndices = node->getColumnIndices(); //копируем индексы

    //bool hasDependencySubset = false;
    //bool hasNonDependencySuperset = false;
    //bool hasUncheckedSubset = false;
    //bool hasUncheckedSuperset = false;

    //bool allSubsetsAreNonDeps = false;

    for (size_t index = 0; index < columnIndices.size(); index++) {
        if (columnIndices[index]) {
            columnIndices[index] = false; //убираем одну из колонок
            auto const subsetNodeIter = this->find(Vertical(node->getSchema(), columnIndices));

            if (subsetNodeIter != this->end() &&
                (subsetNodeIter->second == NodeCategory::minimalDependency ||
                 subsetNodeIter->second == NodeCategory::dependency)
            ) {
                //получили кандидата в минимальную зависимость
                return true;
            }
            columnIndices[index] = true; //возвращаем как было

        } else {
            columnIndices[index] = true;
            auto const supersetNodeIter = this->find(Vertical(node->getSchema(), columnIndices));

            if (supersetNodeIter != this->end() &&
                (supersetNodeIter->second == NodeCategory::maximalNonDependency ||
                 supersetNodeIter->second == NodeCategory::nonDependency)
            ) {
                //получили кандидата в макс. независимость
                return true;
            }
            columnIndices[index] = false; //возвращаем как было
        }
    }

    return false;*/
}

bool LatticeObservations::isCandidate(const shared_ptr<Vertical> &node) {
    auto nodeIter = this->find(*node);
    if (nodeIter == this->end()) {
        return false;
    } else {
        return nodeIter->second == NodeCategory::candidateMaximalNonDependency ||
               nodeIter->second == NodeCategory::candidateMinimalDependency;
    }
}

LatticeObservations::vertical_set LatticeObservations::getUncheckedSubsets(const shared_ptr<Vertical> &node) {
    vertical_set uncheckedSubsets;

    for (auto& subsetNode : node->getParents()) {
        if (this->find(*subsetNode) == this->end()) {
            uncheckedSubsets.insert(std::move(subsetNode));
        }
    }

    return uncheckedSubsets;
}

LatticeObservations::vertical_set LatticeObservations::getUncheckedSupersets(const shared_ptr<Vertical> &node) {
    dynamic_bitset<> indices = node->getColumnIndices();
    vertical_set uncheckedSupersets;

    for (size_t index = 0; index < indices.size(); index++) {
        if (!indices[index]) {
            indices[index] = true;
            auto supersetNode = std::make_shared<Vertical>(node->getSchema(), indices);

            if (this->find(*supersetNode) == this->end()) {
                uncheckedSupersets.insert(std::move(supersetNode));
            }

            indices[index] = false;
        }
    }
    return uncheckedSupersets;
}


