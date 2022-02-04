//
// Created by kek on 16.08.2019.
//

//#include "ColumnCombination.h"

//#include "RelationalSchema.h"

/*ColumnCombination::ColumnCombination(dynamic_bitset<> column_indices_, shared_ptr<RelationalSchema> schema):
    Vertical(){
    this->column_indices_ = std::move(column_indices_);
    this->schema = schema;
}

//full of errors - get rid of it
string ColumnCombination::ToString() {
    string ans = "[";
    string separator;
    auto relation = schema.lock();
    for (unsigned long index = column_indices_.find_first(); index < column_indices_.size(); index = column_indices_.find_next(index)){
        ans += separator + relation->getColumn(index)->GetName();
        separator = ", ";
    }
    return ans;
}*/