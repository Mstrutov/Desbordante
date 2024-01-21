#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "model/table/tuple_index.h"

namespace algos::order {

struct SortedPartition {
    using EquivalenceClasses = std::vector<std::unordered_set<model::TupleIndex>>;
    using PartitionIndex = unsigned long;

    EquivalenceClasses sorted_partition;
    std::unordered_map<unsigned long int, unsigned long int> hash_partition;
    unsigned long num_rows = 0;

    SortedPartition() = default;
    SortedPartition(unsigned long num_rows) : num_rows(num_rows){};
    SortedPartition(EquivalenceClasses&& eq_classes, unsigned long num_rows)
        : sorted_partition(std::move(eq_classes)), num_rows(num_rows){};

    void BuildHashTable();
    void Intersect(SortedPartition const& other);
};

}  // namespace algos::order
