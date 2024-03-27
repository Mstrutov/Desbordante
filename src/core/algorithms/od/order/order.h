#pragma once

#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "algorithms/algorithm.h"
#include "config/tabular_data/input_table_type.h"
#include "dependency_checker.h"
#include "list_lattice.h"
#include "model/table/column_layout_typed_relation_data.h"
#include "order_utility.h"
#include "sorted_partitions.h"

namespace algos::order {

class Order : public Algorithm {
public:
    using SortedPartitions = std::unordered_map<Node, SortedPartition, Hash>;
    using TypedRelation = model::ColumnLayoutTypedRelationData;

    config::InputTable input_table_;
    std::unique_ptr<TypedRelation> typed_relation_;
    SortedPartitions sorted_partitions_;
    std::vector<AttributeList> single_attributes_;
    CandidateSets previous_candidate_sets_;
    CandidateSets candidate_sets_;
    OrderDependencies valid_;
    OrderDependencies merge_invalidated_;
    std::unique_ptr<ListLattice> lattice_;

    void RegisterOptions();
    void LoadDataInternal() override;
    void ResetState() override;
    void PruneSingleEqClassPartitions();
    void CreateSingleColumnSortedPartitions();
    void CreateSortedPartitionsFromSingletons(AttributeList const& attr_list);
    void ComputeDependencies(ListLattice::LatticeLevel const& lattice_level);
    std::vector<AttributeList> Extend(AttributeList const& lhs, AttributeList const& rhs) const;
    bool IsMinimal(AttributeList const& a) const;
    void UpdateCandidateSets();
    void MergePrune();
    void PrintValidOD();
    unsigned long long ExecuteInternal() final;

public:
    OrderDependencies const& GetValidODs() const {
        return valid_;
    }

    Order();
};

}  // namespace algos::order
