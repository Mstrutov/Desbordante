//
// Created by kek on 17.07.19.
//

#include "model/RelationalSchema.h"

#include <utility>

#include "model/ColumnCombination.h"

using namespace std;

RelationalSchema::RelationalSchema(string name, bool isNullEqNull) :
        columns(),
        name(std::move(name)),
        isNullEqNull(isNullEqNull),
        emptyVertical() {
}
// good practice is using std::make_shared instead
shared_ptr<RelationalSchema> RelationalSchema::create(string name, bool isNullEqNull) {
    auto schema = shared_ptr<RelationalSchema>(new RelationalSchema(std::move(name), isNullEqNull));
    schema->init();
    return schema;
}

// this is hard to comprehend
void RelationalSchema::init() {
    emptyVertical.reset(new Vertical(std::move(Vertical::emptyVertical(shared_from_this()))));
}

//TODO: В оригинале тут какая-то срань
Vertical RelationalSchema::getVertical(dynamic_bitset<> indices) {
    if (indices.empty()) return *(this->emptyVertical);

    if (indices.count() == 1){
        return *(std::static_pointer_cast<Vertical>(this->columns[indices.find_first()]));          //TODO: TEMPORAL KOSTYL'
    }

    return ColumnCombination(indices, shared_from_this());
}

string RelationalSchema::getName() { return name; }

vector<shared_ptr<Column>> RelationalSchema::getColumns() { return columns; }

//TODO: assert'ы пофиксить на нормальные эксепшены
shared_ptr<Column> RelationalSchema::getColumn(const string &colName) {
    for (auto &column : columns){
        if (column->name == colName)
            return column;
    }
    assert(0);
}

shared_ptr<Column> RelationalSchema::getColumn(int index) {
    return columns[index];
}

void RelationalSchema::appendColumn(const string& colName) {
    columns.push_back(make_shared<Column>(shared_from_this(), colName, columns.size()));
}

// if you have nothing else to do: push_back through move semantics
void RelationalSchema::appendColumn(shared_ptr<Column> column) {
    columns.push_back(column);
}

int RelationalSchema::getNumColumns() {
    return columns.size();
}

bool RelationalSchema::isNullEqualNull() { return isNullEqNull; }

