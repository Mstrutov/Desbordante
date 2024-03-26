#include "order_utility.h"

#include <iostream>

#include <easylogging++.h>

#include "list_lattice.h"
#include "model/table/column_index.h"

namespace algos::order {

void PrintOD(AttributeList const& lhs, AttributeList const& rhs) {
    for (model::ColumnIndex attr : lhs) {
        LOG(DEBUG) << attr << " ";
    }
    LOG(DEBUG) << "-> ";
    for (model::ColumnIndex attr : rhs) {
        LOG(DEBUG) << attr << " ";
    }
}

std::vector<AttributeList> GetPrefixes(Node const& node) {
    std::vector<AttributeList> res;
    res.reserve(node.size() - 1);
    for (size_t i = 1; i < node.size(); ++i) {
        res.emplace_back(node.begin(), node.begin() + i);
    }
    return res;
}

AttributeList MaxPrefix(AttributeList const& attribute_list) {
    return AttributeList(attribute_list.begin(), attribute_list.end() - 1);
}

bool InUnorderedMap(OrderDependencies const& map, AttributeList const& lhs,
                    AttributeList const& rhs) {
    auto it = map.find(lhs);
    if (it == map.end()) {
        return false;
    }
    if (it->second.find(rhs) == it->second.end()) {
        return false;
    }
    return true;
}

bool AreDisjoint(AttributeList const& a, AttributeList const& b) {
    for (model::ColumnIndex a_atr : a) {
        for (model::ColumnIndex b_atr : b) {
            if (a_atr == b_atr) {
                return false;
            }
        }
    }
    return true;
}

bool StartsWith(AttributeList const& rhs_candidate, AttributeList const& rhs) {
    for (model::ColumnIndex i = 0; i < rhs.size(); ++i) {
        if (rhs[i] != rhs_candidate[i]) {
            return false;
        }
    }
    return true;
}

}  // namespace algos::order
