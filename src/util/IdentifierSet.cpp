#include "IdentifierSet.h"

IdentifierSet::IdentifierSet(ColumnLayoutRelationData const* relation,
                             int index) : relation_(relation), tuple_index_(index) {
    data_.reserve(relation_->getNumColumns());
    for (ColumnData const& col : relation_->getColumnData()) {
        data_.push_back({ .attribute     = col.getColumn(),
                          .cluster_index = col.getProbingTableValue(tuple_index_) });
    }
}

Vertical IdentifierSet::intersect(IdentifierSet const& other) const {
    boost::dynamic_bitset<> intersection(relation_->getNumColumns());
    auto p = data_.begin();
    auto q = other.data_.begin();

    while (p != data_.end() && q != other.data_.end()) {
        if (p->attribute->getIndex() < q->attribute->getIndex()) {
            ++p;
        } else {
            if (q->attribute->getIndex() == p->attribute->getIndex() &&
                p->cluster_index != 0 &&
                p->cluster_index == q->cluster_index) {
                intersection.set(p->attribute->getIndex());
                ++p;
            }
            ++q;
        }
    }

    return relation_->getSchema()->getVertical(intersection);
}

std::string IdentifierSet::toString() const {
    if (data_.empty())
        return "[]";
    std::string str = "[";
    for (auto p = data_.begin(); p != data_.end() - 1; ++p)
        str += "(" + p->attribute->getName() + ", " + std::to_string(p->cluster_index) + "), ";
    str += "(" + data_.back().attribute->getName() + ", " + std::to_string(data_.back().cluster_index) + ")]";
    return str;
}