#include "dependency_checker.h"

#include <algorithm>
#include <unordered_set>

#include "model/table/tuple_index.h"

namespace algos::order {

namespace {
bool SubsetSetDifference(SortedPartition::EquivalenceClass const& a,
                         SortedPartition::EquivalenceClass& b) {
    auto const not_found = b.end();
    for (unsigned long element : a) {
        auto it = b.find(element);
        if (it == not_found) {
            return false;
        }
        b.erase(it);
    }
    return true;
}
}  // namespace

ValidityType CheckForSwap(SortedPartition const& l, SortedPartition const& r) {
    ValidityType res = ValidityType::valid;
    size_t l_i = 0, r_i = 0;
    bool next_l = true, next_r = true;
    SortedPartition::EquivalenceClass l_eq_class;
    SortedPartition::EquivalenceClass r_eq_class;
    SortedPartition::EquivalenceClasses const& l_classes = l.GetEqClasses();
    SortedPartition::EquivalenceClasses const& r_classes = r.GetEqClasses();
    while (l_i < l_classes.size() && r_i < r_classes.size()) {
        if (next_l) {
            l_eq_class = l_classes[l_i];
        }
        if (next_r) {
            r_eq_class = r_classes[r_i];
        }
        if (l_eq_class.size() < r_eq_class.size()) {
            if (!SubsetSetDifference(l_eq_class, r_eq_class)) {
                return ValidityType::swap;
            } else {
                res = ValidityType::merge;
                ++l_i;
                next_l = true;
                next_r = false;
            }
        } else {
            if (!SubsetSetDifference(r_eq_class, l_eq_class)) {
                return ValidityType::swap;
            } else {
                ++r_i;
                next_r = true;
                if (l_eq_class.empty()) {
                    ++l_i;
                    next_l = true;
                } else {
                    next_l = false;
                }
            }
        }
    }
    return res;
}

}  // namespace algos::order
