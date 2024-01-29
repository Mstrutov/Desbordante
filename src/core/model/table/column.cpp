#include "column.h"

#include "model/table/relational_schema.h"
#include "vertical.h"

using namespace std;

bool Column::operator<(Column const& rhs) const {
    // use assert here to check if Columns belong to the same schema?
    return index_ > rhs.index_ && *schema_ == *rhs.schema_;
}

bool Column::operator>(Column const& rhs) const {
    return index_ < rhs.index_ && *schema_ == *rhs.schema_;
}

bool Column::operator==(Column const& rhs) const {
    if (this == &rhs) return true;
    return index_ == rhs.index_ && *schema_ == *rhs.schema_;
}

bool Column::operator!=(Column const& rhs) const {
    return !(*this == rhs);
}

Column::operator Vertical() const {
    return Vertical(*this);
}
