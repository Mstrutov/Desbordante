//
// Created by kek on 18.07.19.
//


#include <utility>

#include "Column.h"
#include "Vertical.h"


using namespace std;

bool Column::operator==(const Column &rhs) {
    if (this == &rhs) return true;
    return index == rhs.index && schema.lock().get() == rhs.schema.lock().get();
}

Column::operator Vertical() const {
    return Vertical(*this);
}
