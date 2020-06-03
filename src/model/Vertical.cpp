//
// Created by kek on 17.07.19.
//

#include "model/Vertical.h"

#include <iostream>
#include <utility>



using namespace std;

Vertical::Vertical(shared_ptr<RelationalSchema> relSchema, dynamic_bitset<> const & indices) :
    columnIndices(indices),
    schema(relSchema) {}

Vertical::Vertical(Column & col) : schema(col.getSchema()){
    columnIndices = dynamic_bitset<>(schema.lock()->getNumColumns());
    columnIndices.set(col.getIndex());
}

Vertical::Vertical(Vertical &&other) noexcept :
    columnIndices(std::move(other.columnIndices)),
    schema(std::move(other.schema)) {
}

Vertical& Vertical::operator=(Vertical &&rhs) noexcept {
    columnIndices = std::move(rhs.columnIndices);
    schema = std::move(rhs.schema);
    return *this;
}

dynamic_bitset<> Vertical::getColumnIndices() const { return columnIndices; }

shared_ptr<RelationalSchema> Vertical::getSchema() { return schema.lock(); }

//TODO: перепроверь цикл
//TODO: перепроверь все операции с Джавой
bool Vertical::contains(Vertical &that) {
    dynamic_bitset<>& thisIndices = columnIndices;
    dynamic_bitset<>& thatIndices = that.columnIndices;
    if(thisIndices.size() < thatIndices.size()) return false;
    for (unsigned long columnIndex = thatIndices.find_first(); columnIndex < thatIndices.size(); columnIndex = thatIndices.find_next(columnIndex)){
        if (!(thisIndices[columnIndex])) return false;
    }
    return true;
}

bool Vertical::intersects(Vertical &that) {
    dynamic_bitset<>& thisIndices = columnIndices;
    dynamic_bitset<>& thatIndices = that.columnIndices;
    return thisIndices.intersects(thatIndices);
}

Vertical Vertical::Union(Vertical const &that) {
    dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices |= that.columnIndices;
    return schema.lock()->getVertical(retainedColumnIndices);
}

Vertical Vertical::project(Vertical &that) {
    dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices &= that.columnIndices;
    return schema.lock()->getVertical(retainedColumnIndices);
}

//TODO: check
Vertical Vertical::without(Vertical const & that) const {
    dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices &= ~that.columnIndices;
    //retainedColumnIndices = ~retainedColumnIndices;
    return schema.lock()->getVertical(retainedColumnIndices);
}

//TODO: UNUSED METHOD - CHECK ITS VALIDITY
Vertical Vertical::invert() {
    auto relation = schema.lock();
    dynamic_bitset<> flippedIndices(columnIndices);
    flippedIndices.resize(relation->getNumColumns());
    flippedIndices.flip();
    return relation->getVertical(flippedIndices);
}

Vertical Vertical::invert(Vertical &scope) {
    dynamic_bitset<> flippedIndices(columnIndices);
    flippedIndices ^= scope.columnIndices;
    return schema.lock()->getVertical(flippedIndices);
}

Vertical Vertical::emptyVertical(shared_ptr<RelationalSchema> relSchema) {
    return Vertical(relSchema, dynamic_bitset<>(relSchema->getNumColumns()));
}

int Vertical::getArity() const {
    return columnIndices.count();
}

vector<shared_ptr<Column>> Vertical::getColumns() const {
    //dynamic_bitset<> returnColumnIndices = getColumnIndices();
    auto relation = schema.lock();
    vector<shared_ptr<Column>> columns;
    for (int index = columnIndices.find_first();
         index != -1;//dynamic_bitset<>::npos;
         index = columnIndices.find_next(index)) {
        columns.push_back(relation->getColumns()[index]);
    }
    return columns;
}

string Vertical::toString() const {
    string result = "[";

    if ((int)columnIndices.find_first() == -1)
        return "Empty Vertical";

    auto relation = schema.lock();
    int i = columnIndices.find_next(0);
    for (int index = columnIndices.find_first();
         index != -1;//dynamic_bitset<>::npos;
         index = columnIndices.find_next(index)) {
        result += relation->getColumn(index)->getName() + " ";   //to_string(index) + " ";
    }
    result[result.size() - 1] = ']';
    return result;
}

vector<shared_ptr<Vertical>> Vertical::getParents() {
    if (getArity() < 2) return vector<shared_ptr<Vertical>>();
    vector<shared_ptr<Vertical>> parents(getArity());
    int i = 0;
    for (size_t columnIndex = columnIndices.find_first();
         columnIndex != dynamic_bitset<>::npos;
         columnIndex = columnIndices.find_next(columnIndex)) {
        columnIndices.reset(columnIndex);
        parents[i++] = std::make_shared<Vertical>(getSchema()->getVertical(columnIndices));
        columnIndices.set(columnIndex);
    }
    return parents;
}
