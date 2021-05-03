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
#include "CustomComparator.h"

class DFD : public FDAlgorithm {
private:
    LatticeObservations observations;
    shared_ptr<PartitionStorage> partitionStorage;
    DependenciesMap dependenciesMap;
    NonDependenciesMap nonDependenciesMap;
    shared_ptr<ColumnLayoutRelationData> relation;

    std::unordered_set<shared_ptr<Vertical>> minimalDeps; //TODO мб их определять либо в функции execute, либо полями класса
    std::unordered_set<shared_ptr<Vertical>> maximalNonDeps;

    std::stack<shared_ptr<Vertical>> trace;

    std::random_device rd; // для генерации случайных чисел
    std::mt19937 gen;

    void findLHSs(shared_ptr<Column const> const& rhs, shared_ptr<RelationalSchema> schema); //TODO: нужен ли второй параметр?; мб переименовать типа findDeps
    shared_ptr<Vertical> pickNextNode(shared_ptr<Vertical> const& node);
    std::list<shared_ptr<Vertical>> generateNextSeeds(shared_ptr<Column const> const& currentRHS);
    shared_ptr<Vertical> takeRandom(std::list<shared_ptr<Vertical>> & nodeList);
    shared_ptr<Vertical> takeRandom(std::vector<shared_ptr<Vertical>> const& nodeList);
    void minimize(std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator> & nodeList);

public:
    explicit DFD(std::filesystem::path const& path, char separator = ',', bool hasHeader = true);

    unsigned long long execute() override;
};
