//
// Created by alex on 13.04.2021.
//

#include "LatticeObservations.h"

bool LatticeObservations::checkIfMinimalDependency(shared_ptr<Vertical> const& vertical) {
    //бежим по множеству подмножеств и смотрим. если получили все независимости, то это мин зависимость
    //TODO что если узел единичный
    dynamic_bitset<> columnIndices = vertical->getColumnIndices(); //копируем
    for (size_t index = columnIndices.find_first(); index < columnIndices.size(); index = columnIndices.find_next(index)) {
        columnIndices[index] = false; //убираем одну из колонок
        auto subsetVertical = this->find(Vertical(vertical->getSchema(), columnIndices));

        //если какое-то подмножество не посещено либо оно не является антизависимостью, то не подходит
        //TODO может быть такое говно что если нашли какую-то неантизависимость, то убираем кандидатство?
        if (subsetVertical == this->end() ||
           (subsetVertical->second != NodeCategory::nonDependency &&
            subsetVertical->second != NodeCategory::maximalNonDependency)
        ) {
            return false;
        }

        //если все норм
        columnIndices[index] = true; //возвращаем как было
    }
    return true;
}


