//
// Created by alex on 13.04.2021.
//

#include "LatticeObservations.h"

NodeCategory LatticeObservations::updateDependencyCategory(const shared_ptr<Vertical> &vertical) {
    //бежим по множеству подмножеств и смотрим. если получили все независимости, то это мин зависимость
    //TODO что если узел единичный
    dynamic_bitset<> columnIndices = vertical->getColumnIndices(); //копируем индексы
    for (size_t index = columnIndices.find_first(); index < columnIndices.size(); index = columnIndices.find_next(index)) {
        columnIndices[index] = false; //убираем одну из колонок
        auto const subsetVerticalIter = this->find(Vertical(vertical->getSchema(), columnIndices)); //TODO передаем временнй объект??

        //если какое-то подмножество не посещено либо оно не является антизависимостью, то не подходит
        //TODO может быть такое говно что если нашли какую-то неантизависимость, то убираем кандидатство?

        if (subsetVerticalIter == this->end()) {
            //если нашли нерассмотренное подмножество
            return NodeCategory::candidateMinimalDependency;
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
    return NodeCategory::minimalDependency;
}

NodeCategory LatticeObservations::updateNonDependencyCategory(const shared_ptr<Vertical> &vertical) {
    //копируем индексы колонок, которые не содержатся в искомом вертикале, чтобы исследовать надмножетва
    dynamic_bitset<> invertedColumnIndices = vertical->getColumnIndices().flip();

    for (size_t index = invertedColumnIndices.find_first(); index < invertedColumnIndices.size(); index = invertedColumnIndices.find_next(index)) {
        invertedColumnIndices[index] = false; //убираем одну из колонок
        auto const supersetVerticalIter = this->find(Vertical(vertical->getSchema(), invertedColumnIndices.flip())); //TODO !!!лучше переделать без второго flip а просто бежать циклом по нулям

        if (supersetVerticalIter == this->end()) {
            //если нашли нерассмотренное надмножество
            return NodeCategory::candidateMaximalNonDependency;
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
    return NodeCategory::maximalNonDependency;
}


