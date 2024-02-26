#pragma once

#include <list>
#include <random>
#include <stack>
#include <vector>

#include "../column_order/column_order.h"
#include "../lattice_observations/lattice_observations.h"
#include "../partition_storage/partition_storage.h"
#include "../pruning_maps/dependencies_map.h"
#include "../pruning_maps/non_dependencies_map.h"
#include "model/table/vertical.h"

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
