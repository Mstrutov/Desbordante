#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "model/table/tuple_index.h"

namespace algos::order {

class SortedPartition {
public:
    using EquivalenceClass = std::unordered_set<model::TupleIndex>;
    using EquivalenceClasses = std::vector<EquivalenceClass>;
    using PartitionIndex = unsigned long;
    using HashProduct = std::unordered_map<PartitionIndex, SortedPartition::EquivalenceClasses>;

private:
    EquivalenceClasses sorted_partition_;
    std::unordered_map<unsigned long int, unsigned long int> hash_partition_;
    unsigned long num_rows_ = 0;

    void BuildHashTable();
    HashProduct BuildHashProduct(SortedPartition const& other);

public:
    SortedPartition() = default;
    SortedPartition(unsigned long num_rows) : num_rows_(num_rows){};
    SortedPartition(EquivalenceClasses&& eq_classes, unsigned long num_rows)
        : sorted_partition_(std::move(eq_classes)), num_rows_(num_rows){};
    void Intersect(SortedPartition const& other);

    EquivalenceClasses const& GetEqClasses() const {
        return sorted_partition_;
    }

    std::size_t Size() const {
        return sorted_partition_.size();
    }
};

}  // namespace algos::order
