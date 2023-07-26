#pragma once

#include <unordered_set>
#include <vector>

namespace algos::order {

struct SortedPartition {
    using EquivalenceClasses = std::vector<std::unordered_set<unsigned long>>;

    EquivalenceClasses sorted_partition;

    SortedPartition() = default;
    SortedPartition(EquivalenceClasses&& eq_classes) : sorted_partition(std::move(eq_classes)){};
};

}  // namespace algos::order
