#pragma once

#include <memory>
#include <unordered_map>

#include <boost/container_hash/hash.hpp>

#include "algorithms/algorithm.h"
#include "model/table/column_layout_typed_relation_data.h"
#include "config/tabular_data/input_table_type.h"
#include "sorted_partitions.h"

namespace algos::order {

class Order : public Algorithm {
private:
    using TypedRelation = model::ColumnLayoutTypedRelationData;
    using SortedPartitions = std::unordered_map<std::vector<unsigned int>, SortedPartition,
                                                boost::hash<std::vector<unsigned int>>>;

    config::InputTable input_table_;
    std::unique_ptr<TypedRelation> typed_relation_;
    SortedPartitions sorted_partitions_;

    void RegisterOptions();
    void LoadDataInternal() override;
    void ResetState() override;
    void CreateSortedPartitions();
    unsigned long long ExecuteInternal() final;

public:
    Order();
};

}  // namespace algos::order
