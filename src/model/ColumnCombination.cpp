//
// Created by kek on 16.08.2019.
//

#include "model/ColumnCombination.h"

#include "model/RelationalSchema.h"

ColumnCombination::ColumnCombination(dynamic_bitset<> columnIndices, shared_ptr<RelationalSchema> schema):
    Vertical(){
    this->columnIndices = std::move(columnIndices);
    this->schema = schema;
}

//full of errors - get rid of it
string ColumnCombination::toString() {
    string ans = "[";
    string separator;
    auto relation = schema.lock();
    for (unsigned long index = columnIndices.find_first(); index < columnIndices.size(); index = columnIndices.find_next(index)){
        ans += separator + relation->getColumn(index)->getName();
        separator = ", ";
    }
    return ans;
}