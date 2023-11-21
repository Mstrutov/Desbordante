#include "dependency_checker.h"

#include <algorithm>
#include <unordered_set>

namespace algos::order {

namespace {
bool SubsetSetDifference(std::unordered_set<unsigned long> const& a,
                         std::unordered_set<unsigned long>& b) {
    auto const not_found = b.end();
    for (unsigned long element : a) {
        if (b.find(element) == not_found) {
            return false;
        } else {
            b.erase(element);
        }
    }
    return true;
}
}  // namespace

ValidityType CheckForSwap(SortedPartition const& l, SortedPartition const& r) {
    ValidityType res = ValidityType::valid;
    size_t l_i = 0, r_i = 0;
    bool next_l = true, next_r = true;
    std::unordered_set<unsigned long> l_eq_class;
    std::unordered_set<unsigned long> r_eq_class;
    while (l_i < l.sorted_partition.size() && r_i < r.sorted_partition.size()) {
        if (next_l) {
            l_eq_class = l.sorted_partition[l_i];
        }
        if (next_r) {
            r_eq_class = r.sorted_partition[r_i];
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
