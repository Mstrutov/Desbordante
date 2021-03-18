#include <iostream>
#include <utility>

#include "Vertical.h"

Vertical::Vertical(RelationalSchema const* relSchema, boost::dynamic_bitset<> indices) :
    columnIndices(std::move(indices)),
    schema(relSchema) {}


Vertical::Vertical(Column const& col) : schema(col.getSchema()){
    columnIndices = boost::dynamic_bitset<>(schema->getNumColumns());
    columnIndices.set(col.getIndex());
}

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

Vertical Vertical::Union(Vertical const &that) const {
    boost::dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices |= that.columnIndices;
    return schema->getVertical(retainedColumnIndices);
}

Vertical Vertical::Union(Column const& that) const {
    boost::dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices.set(that.getIndex());
    return schema->getVertical(retainedColumnIndices);
}

Vertical Vertical::project(Vertical const& that) const {
    boost::dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices &= that.columnIndices;
    return schema->getVertical(retainedColumnIndices);
}


Vertical Vertical::without(Vertical const& that) const {
    boost::dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices &= ~that.columnIndices;
    return schema->getVertical(retainedColumnIndices);
}

Vertical Vertical::without(Column const& that) const {
    boost::dynamic_bitset<> retainedColumnIndices(columnIndices);
    retainedColumnIndices.reset(that.getIndex());
    return schema->getVertical(retainedColumnIndices);
}

Vertical Vertical::invert() const {
    boost::dynamic_bitset<> flippedIndices(columnIndices);
    flippedIndices.resize(schema->getNumColumns());
    flippedIndices.flip();
    return schema->getVertical(flippedIndices);
}

Vertical Vertical::invert(Vertical const& scope) const {
    boost::dynamic_bitset<> flippedIndices(columnIndices);
    flippedIndices ^= scope.columnIndices;
    return schema->getVertical(flippedIndices);
}

std::unique_ptr<Vertical> Vertical::emptyVertical(RelationalSchema const* relSchema) {
    return std::make_unique<Vertical>(relSchema, boost::dynamic_bitset<>(relSchema->getNumColumns()));
}

std::vector<Column const*> Vertical::getColumns() const {
    std::vector<Column const*> columns;
    for (size_t index = columnIndices.find_first();
         index != boost::dynamic_bitset<>::npos;
         index = columnIndices.find_next(index)) {
        columns.push_back(schema->getColumns()[index].get());
    }
    return columns;
}

std::string Vertical::toString() const {
    std::string result = "[";

    if (columnIndices.find_first() == boost::dynamic_bitset<>::npos)
        return "[]";

    for (size_t index = columnIndices.find_first();
         index != boost::dynamic_bitset<>::npos;
         index = columnIndices.find_next(index)) {
        result += schema->getColumn(index)->getName();
        if (columnIndices.find_next(index) != boost::dynamic_bitset<>::npos) {
            result += ' ';
        }
    }

    result += ']';

    return result;
}

std::string Vertical::toIndicesString() const {
    std::string result = "[";

    if (columnIndices.find_first() == boost::dynamic_bitset<>::npos) {
        return "[]";
    }

    for (size_t index = columnIndices.find_first();
         index != boost::dynamic_bitset<>::npos;
         index = columnIndices.find_next(index)) {
        result += std::to_string(index);
        if (columnIndices.find_next(index) != boost::dynamic_bitset<>::npos) {
            result += ',';
        }
    }

    result += ']';

    return result;
}



std::vector<Vertical> Vertical::getParents() const {
    if (getArity() < 2) return std::vector<Vertical>();
    std::vector<Vertical> parents(getArity());
    int i = 0;
    for (size_t columnIndex = columnIndices.find_first();
         columnIndex != boost::dynamic_bitset<>::npos;
         columnIndex = columnIndices.find_next(columnIndex)) {
        auto parentColumnIndices = columnIndices;
        parentColumnIndices.reset(columnIndex);
        parents[i++] = getSchema()->getVertical(std::move(parentColumnIndices));
    }
    return parents;
}
