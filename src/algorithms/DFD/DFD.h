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

class DFD : public FDAlgorithm {
private:
    //using vertical_set = std::unordered_set<shared_ptr<Vertical>, std::hash<shared_ptr<Vertical>>, custom_comparator>;

    LatticeObservations observations;
    std::unique_ptr<PartitionStorage> partitionStorage;
    DependenciesMap dependenciesMap;
    NonDependenciesMap nonDependenciesMap;
    std::unique_ptr<ColumnLayoutRelationData> relation;

    std::unordered_set<Vertical> minimalDeps; //TODO мб их определять либо в функции execute, либо полями класса
    std::unordered_set<Vertical> maximalNonDeps;

    std::stack<Vertical> trace;

    std::random_device rd; // для генерации случайных чисел
    std::mt19937 gen;

    void findLHSs(Column const* const rhs, RelationalSchema const* const schema); //TODO: нужен ли второй параметр?; мб переименовать типа findDeps
    Vertical pickNextNode(Vertical const &node, size_t rhsIndex);
    std::list<Vertical> generateNextSeeds(Column const* const currentRHS);
    Vertical takeRandom(std::list<Vertical> & nodeList);
    Vertical takeRandom(std::vector<Vertical> const& nodeList);
    Vertical takeRandom(std::unordered_set<Vertical>&);
    void minimize(std::unordered_set<Vertical> &);
    static void substractSets(std::unordered_set<Vertical> & set, std::unordered_set<Vertical> const& setToSubstract);
    bool inferCategory(Vertical const& node, size_t rhsIndex);

public:
    explicit DFD(std::filesystem::path const& path, char separator = ',', bool hasHeader = true);

    unsigned long long execute() override;
};
