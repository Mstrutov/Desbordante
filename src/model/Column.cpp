//
// Created by kek on 18.07.19.
//


#include <utility>

#include "Column.h"
#include "Vertical.h"


using namespace std;


unsigned int Column::getIndex() const { return index;}

string Column::getName() const { return name; }



string Column::toString() const { return "[" + name + "]";}

bool Column::operator==(const Column &rhs) {
    if (this == &rhs) return true;
    return index == rhs.index && schema.lock().get() == rhs.schema.lock().get();
}

std::shared_ptr<RelationalSchema> Column::getSchema() const {
    return schema.lock();
}

Column::operator Vertical() const {
    return Vertical(*this);
}
