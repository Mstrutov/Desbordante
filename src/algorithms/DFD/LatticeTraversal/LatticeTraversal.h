#pragma once

#include <stack>

#include "FDAlgorithm.h"
#include "Vertical.h"
#include "DFD/ColumnOrder/ColumnOrder.h"
#include "DFD/LatticeObservations/LatticeObservations.h"
#include "DFD/PruningMaps/DependenciesMap.h"
#include "DFD/PruningMaps/NonDependenciesMap.h"
#include "DFD/PartitionStorage/PartitionStorage.h"

class LatticeTraversal {
private:
    Column const* const rhs;

    std::unordered_set<Vertical> minimalDeps;
    std::unordered_set<Vertical> maximalNonDeps;
    DependenciesMap dependenciesMap;
    NonDependenciesMap nonDependenciesMap;
    LatticeObservations observations;
    std::stack<Vertical> trace;
    ColumnOrder const columnOrder;

    std::vector<Vertical> const& uniqueColumns;
    ColumnLayoutRelationData const* const relation;
    PartitionStorage * const partitionStorage;

    std::random_device rd;
    std::mt19937 gen;

    bool inferCategory(Vertical const& node, unsigned int rhsIndex);
    Vertical pickNextNode(Vertical const &node, unsigned int rhsIndex);
    std::stack<Vertical> generateNextSeeds(Column const* const currentRHS);

    std::list<Vertical> minimize(std::unordered_set<Vertical> const&) const;
    Vertical const& takeRandom(std::unordered_set<Vertical> &nodeSet);
    static void substractSets(std::unordered_set<Vertical> & set, std::unordered_set<Vertical> const& setToSubstract);

public:
    LatticeTraversal(Column const* const rhs,
                     ColumnLayoutRelationData const* const relation,
                     std::vector<Vertical> const& uniqueVerticals,
                     PartitionStorage * const partitionStorage);

    std::unordered_set<Vertical> findLHSs();
};
