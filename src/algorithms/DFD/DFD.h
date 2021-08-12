//
// Created by alexandrsmirn
//

#pragma once

#include <random>
#include <stack>

#include "FDAlgorithm.h"
#include "DFD/LatticeObservations/LatticeObservations.h"
#include "PLICache.h"
#include "DependenciesMap.h"
#include "NonDependenciesMap.h"
#include "PartitionStorage/PartitionStorage.h"
#include "Vertical.h"
#include "ColumnOrder/ColumnOrder.h"

class DFD : public FDAlgorithm {
private:
    std::unique_ptr<PartitionStorage> partitionStorage;
    std::unique_ptr<ColumnLayoutRelationData> relation;
    DependenciesMap dependenciesMap;
    NonDependenciesMap nonDependenciesMap;
    LatticeObservations observations;
    ColumnOrder columnOrder;

    std::unordered_set<Vertical> minimalDeps;
    std::unordered_set<Vertical> maximalNonDeps;

    std::stack<Vertical> trace;

    std::random_device rd; // для генерации случайных чисел
    std::mt19937 gen;

    void findLHSs(Column const* const rhs);
    bool inferCategory(Vertical const& node, size_t rhsIndex);
    Vertical pickNextNode(Vertical const &node, size_t rhsIndex);
    std::stack<Vertical> generateNextSeeds(Column const* const currentRHS);

    std::list<Vertical> minimize(std::unordered_set<Vertical> const&);
    const Vertical & takeRandom(std::unordered_set<Vertical> &nodeSet);
    static void substractSets(std::unordered_set<Vertical> & set, std::unordered_set<Vertical> const& setToSubstract);

public:
    explicit DFD(std::filesystem::path const& path, char separator = ',', bool hasHeader = true);

    unsigned long long execute() override;
};
