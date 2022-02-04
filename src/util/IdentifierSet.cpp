#include "IdentifierSet.h"

namespace util {

IdentifierSet::IdentifierSet(ColumnLayoutRelationData const* const relation,
                             int index) : relation_(relation), tuple_index_(index) {
    data_.reserve(relation_->GetNumColumns());
    for (ColumnData const& col : relation_->GetColumnData()) {
        data_.push_back({ .attribute     = col.GetColumn(),
                          .cluster_index = col.GetProbingTableValue(tuple_index_) });
    }
}

std::string IdentifierSet::ToString() const {
    if (data_.empty()) {
        return "[]";
    }

    std::string str = "[";
    for (auto p = data_.begin(); p != data_.end() - 1; ++p) {
        str += "(" + p->attribute->GetName() + ", " + std::to_string(p->cluster_index) + "), ";
    }
    str += "(" + data_.back().attribute->GetName() + ", " + std::to_string(data_.back().cluster_index) + ")]";
    return str;
}

} // namespace util

