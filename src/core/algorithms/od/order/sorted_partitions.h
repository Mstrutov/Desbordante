#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace algos::order {

struct SortedPartition {
    using EquivalenceClasses = std::vector<std::unordered_set<unsigned long>>;

    EquivalenceClasses sorted_partition;
    std::unordered_map<unsigned long, unsigned long> hash_partition;

    SortedPartition() = default;
    SortedPartition(EquivalenceClasses&& eq_classes) : sorted_partition(std::move(eq_classes)){};

    void BuildHashTable();
    SortedPartition operator*(SortedPartition const& other);
};

}  // namespace algos::order
