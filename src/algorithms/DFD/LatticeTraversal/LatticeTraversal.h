#pragma once

#include <stack>

#include "Vertical.h"
#include "DFD/ColumnOrder/ColumnOrder.h"
#include "DFD/LatticeObservations/LatticeObservations.h"
#include "DFD/PruningMaps/DependenciesMap.h"
#include "DFD/PruningMaps/NonDependenciesMap.h"
#include "DFD/PartitionStorage/PartitionStorage.h"

class LatticeTraversal {
private:
    Column const* const rhs_;

    std::unordered_set<Vertical> minimal_deps_;
    std::unordered_set<Vertical> maximal_non_deps_;
    DependenciesMap dependencies_map_;
    NonDependenciesMap non_dependencies_map_;
    LatticeObservations observations_;
    std::stack<Vertical> trace_;
    ColumnOrder const column_order_;

    std::vector<Vertical> const& unique_columns_;
    ColumnLayoutRelationData const* const relation_;
    PartitionStorage* const partition_storage_;

    std::random_device rd_;
    std::mt19937 gen_;

    bool InferCategory(Vertical const& node, unsigned int rhs_index);
    Vertical PickNextNode(Vertical const& node, unsigned int rhs_index);
    std::stack<Vertical> GenerateNextSeeds(Column const* const current_rhs);

    std::list<Vertical> Minimize(std::unordered_set<Vertical> const& node_list) const;
    Vertical const& TakeRandom(std::unordered_set<Vertical>& node_set);
    static void SubstractSets(std::unordered_set<Vertical>& set,
                              std::unordered_set<Vertical> const& set_to_substract);

public:
    LatticeTraversal(Column const* const rhs, ColumnLayoutRelationData const* const relation,
                     std::vector<Vertical> const& unique_verticals,
                     PartitionStorage* const partition_storage);

    std::unordered_set<Vertical> FindLHSs();
};
