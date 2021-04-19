//
// Created by kek on 17.07.19.
//

#include "Vertical.h"

#include <iostream>
#include <utility>



using namespace std;

Vertical::Vertical(shared_ptr<RelationalSchema> relSchema, dynamic_bitset<> const & indices) :
    columnIndices(indices),
    schema(relSchema) {}

Vertical::Vertical(Column const& col) : schema(col.getSchema()){
    columnIndices = dynamic_bitset<>(schema.lock()->getNumColumns());
    columnIndices.set(col.getIndex());
}


//TODO: перепроверь цикл
//TODO: перепроверь все операции с Джавой
bool Vertical::contains(Vertical const& that) const {
    boost::dynamic_bitset<> const& thatIndices = that.columnIndices;
    if (columnIndices.size() < thatIndices.size()) return false;

    for (size_t columnIndex = thatIndices.find_first();
         columnIndex < thatIndices.size();
         columnIndex = thatIndices.find_next(columnIndex)) {
        if (!(columnIndices[columnIndex])) return false;
    }
    return true;
}

bool Vertical::contains(Column const& that) const {
    return columnIndices.test(that.getIndex());
}

bool Vertical::intersects(Vertical const& that) const {
    boost::dynamic_bitset<> const& thatIndices = that.columnIndices;
    return columnIndices.intersects(thatIndices);
}

std::shared_ptr<Vertical> Vertical::Union(Vertical const& that) const {
    boost::dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices |= that.columnIndices;
    return schema.lock()->getVertical(retainedColumnIndices);
}

std::shared_ptr<Vertical> Vertical::Union(Column const& that) const {
    boost::dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices.set(that.getIndex());
    return schema.lock()->getVertical(retainedColumnIndices);
}

std::shared_ptr<Vertical> Vertical::project(Vertical const& that) const {
    boost::dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices &= that.columnIndices;
    return schema.lock()->getVertical(retainedColumnIndices);
}

std::shared_ptr<Vertical> Vertical::without(Vertical const& that) const {
    dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices &= ~that.columnIndices;
    return schema.lock()->getVertical(retainedColumnIndices);
}

std::shared_ptr<Vertical> Vertical::without(Column const& that) const {
    dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices.reset(that.getIndex());
    return schema.lock()->getVertical(retainedColumnIndices);
}

std::shared_ptr<Vertical> Vertical::invert() const {
    auto relation = schema.lock();
    dynamic_bitset<> flippedIndices(columnIndices);
    flippedIndices.resize(relation->getNumColumns());
    flippedIndices.flip();
    return relation->getVertical(flippedIndices);
}

std::shared_ptr<Vertical> Vertical::invert(Vertical const& scope) const {
    dynamic_bitset<> flippedIndices(columnIndices);
    flippedIndices ^= scope.columnIndices;
    return schema.lock()->getVertical(flippedIndices);
}

Vertical Vertical::emptyVertical(shared_ptr<RelationalSchema> relSchema) {
    return Vertical(relSchema, dynamic_bitset<>(relSchema->getNumColumns()));
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

    if (columnIndices.find_first() == boost::dynamic_bitset<>::npos)
        return "[]";

    auto relation = schema.lock();
    for (size_t index = columnIndices.find_first();
         index != dynamic_bitset<>::npos;
         index = columnIndices.find_next(index)) {
        result += relation->getColumn(index)->getName();
        if (columnIndices.find_next(index) != dynamic_bitset<>::npos) {
            result += ' ';
        }
    }

    result += ']';

    return result;
}

std::string Vertical::toIndicesString() const {
    string result = "[";

    if (columnIndices.find_first() == boost::dynamic_bitset<>::npos) {
        return "[]";
    }

    auto relation = schema.lock();
    for (size_t index = columnIndices.find_first();
         index != dynamic_bitset<>::npos;
         index = columnIndices.find_next(index)) {
        result += std::to_string(index);
        if (columnIndices.find_next(index) != dynamic_bitset<>::npos) {
            result += ',';
        }
    }

    result += ']';

    return result;
}


vector<shared_ptr<Vertical>> Vertical::getParents() {
    if (getArity() < 2) return vector<shared_ptr<Vertical>>();
    vector<shared_ptr<Vertical>> parents(getArity());
    int i = 0;
    for (size_t columnIndex = columnIndices.find_first();
         columnIndex != dynamic_bitset<>::npos;
         columnIndex = columnIndices.find_next(columnIndex)) {
        auto parentColumnIndices = columnIndices;
        parentColumnIndices.reset(columnIndex);
        parents[i++] = getSchema()->getVertical(std::move(parentColumnIndices));
    }
    return parents;
}
