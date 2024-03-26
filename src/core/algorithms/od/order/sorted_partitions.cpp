#include "sorted_partitions.h"

#include <unordered_map>
#include <vector>

#include "model/table/tuple_index.h"

namespace algos::order {

void SortedPartition::BuildHashTable() {
    hash_partition_.reserve(num_rows_);
    for (PartitionIndex i = 0; i < sorted_partition_.size(); ++i) {
        if (sorted_partition_[i].size() == 1) {
            continue;
        }
        for (model::TupleIndex tuple_index : sorted_partition_[i]) {
            hash_partition_.emplace(tuple_index, i);
        }
    }
}

SortedPartition::HashProduct SortedPartition::BuildHashProduct(SortedPartition const& other) {
    HashProduct hash_product;
    hash_product.reserve(hash_partition_.size());
    for (EquivalenceClass const& eq_class : other.sorted_partition_) {
        if (other.sorted_partition_.size() <= 1) {
            break;
        }
        std::unordered_set<PartitionIndex> visited_positions;
        for (model::TupleIndex tuple_index : eq_class) {
            auto it = hash_partition_.find(tuple_index);
            if (it == hash_partition_.end()) {
                continue;
            }
            PartitionIndex position = it->second;
            visited_positions.insert(position);
            auto pos_it = hash_product.find(position);
            if (pos_it == hash_product.end()) {
                pos_it = hash_product.emplace(position, EquivalenceClasses{{}}).first;
            }
            pos_it->second.back().insert(tuple_index);
        }
        for (PartitionIndex position : visited_positions) {
            hash_product[position].push_back({});
        }
    }
    return hash_product;
}

void SortedPartition::Intersect(SortedPartition const& other) {
    BuildHashTable();
    HashProduct hash_product = BuildHashProduct(other);
    SortedPartition res(num_rows_);
    res.sorted_partition_.reserve(num_rows_);
    for (std::size_t i = 0; i < sorted_partition_.size(); ++i) {
        if (EquivalenceClass const& eq_class = sorted_partition_[i]; eq_class.size() == 1) {
            res.sorted_partition_.push_back(eq_class);
        } else {
            for (EquivalenceClass const& eq_class : hash_product[i]) {
                if (!eq_class.empty()) {
                    res.sorted_partition_.push_back(std::move(eq_class));
                }
            }
        }
    }
    res.sorted_partition_.shrink_to_fit();
    *this = std::move(res);
}

}  // namespace algos::order
