#include "IdentifierSet.h"

IdentifierSet::IdentifierSet(ColumnLayoutRelationData const* const relation,
                             int index) : relation_(relation), tuple_index_(index) {
    data_.reserve(relation_->getNumColumns());
    for (ColumnData const& col : relation_->getColumnData()) {
        data_.push_back({ .attribute     = col.getColumn(),
                          .cluster_index = col.getProbingTableValue(tuple_index_) });
    }
}

std::string IdentifierSet::toString() const {
    if (data_.empty()) {
        return "[]";
    }

    std::string str = "[";
    for (auto p = data_.begin(); p != data_.end() - 1; ++p) {
        str += "(" + p->attribute->getName() + ", " + std::to_string(p->cluster_index) + "), ";
    }
    str += "(" + data_.back().attribute->getName() + ", " + std::to_string(data_.back().cluster_index) + ")]";
    return str;
}
