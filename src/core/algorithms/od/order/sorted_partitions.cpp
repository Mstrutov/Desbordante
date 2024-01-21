#include "sorted_partitions.h"

#include <unordered_map>
#include <vector>

#include "model/table/tuple_index.h"

namespace algos::order {

void SortedPartition::BuildHashTable() {
    hash_partition.reserve(num_rows);
    for (PartitionIndex i = 0; i < sorted_partition.size(); ++i) {
        if (sorted_partition[i].size() == 1) {
            continue;
        }
        for (model::TupleIndex tuple_index : sorted_partition[i]) {
            hash_partition[tuple_index] = i;
        }
    }
}

void SortedPartition::Intersect(SortedPartition const& other) {
    BuildHashTable();
    std::unordered_map<PartitionIndex, SortedPartition::EquivalenceClasses> hash_product;
    hash_product.reserve(hash_partition.size());
    for (std::unordered_set<model::TupleIndex> const& eq_class : other.sorted_partition) {
        if (other.sorted_partition.size() <= 1) {
            break;
        }
        std::unordered_set<PartitionIndex> visited_positions;
        for (model::TupleIndex tuple_index : eq_class) {
            if (hash_partition.find(tuple_index) == hash_partition.end()) {
                continue;
            }
            PartitionIndex position = hash_partition[tuple_index];
            visited_positions.insert(position);
            if (hash_product.find(position) == hash_product.end()) {
                hash_product[position] = {{}};
            }
            hash_product[position].back().insert(tuple_index);
        }
        for (PartitionIndex position : visited_positions) {
            hash_product[position].push_back({});
        }
    }
    SortedPartition res(num_rows);
    res.sorted_partition.reserve(num_rows);
    for (std::size_t i = 0; i < sorted_partition.size(); ++i) {
        if (sorted_partition[i].size() == 1) {
            res.sorted_partition.push_back(sorted_partition[i]);
        } else {
            for (std::unordered_set<model::TupleIndex> const& eq_class : hash_product[i]) {
                if (!eq_class.empty()) {
                    res.sorted_partition.push_back(std::move(eq_class));
                }
            }
        }
    }
    res.sorted_partition.shrink_to_fit();
    *this = std::move(res);
}

}  // namespace algos::order
