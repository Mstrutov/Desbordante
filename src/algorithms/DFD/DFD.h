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

class DFD : public FDAlgorithm {
private:
    LatticeObservations observations;
    shared_ptr<PLICache> partitionCache;
    DependenciesMap dependencies;
    NonDependenciesMap nonDependenciesMap;
    shared_ptr<ColumnLayoutRelationData> relation;

    std::vector<shared_ptr<Vertical>> minimalDeps; //TODO мб их определять либо в функции execute, либо полями класса
    std::vector<shared_ptr<Vertical>> maximalNonDeps;

    std::stack<shared_ptr<Vertical>> trace;

    std::random_device rd; // для генерации случайных чисел
    std::mt19937 gen;

    void findLHSs(shared_ptr<Column> rhs, shared_ptr<RelationalSchema> schema); //TODO: нужен ли второй параметр?; мб переименовать типа findDeps
    shared_ptr<Vertical> pickNextNode(shared_ptr<Vertical> node);
    std::list<shared_ptr<Vertical>> generateNextSeeds();
    shared_ptr<Vertical> takeRandom(std::list<shared_ptr<Vertical>> const& nodeList);
    shared_ptr<Vertical> takeRandom(std::vector<shared_ptr<Vertical>> const& nodeList);
    void minimize(std::unordered_set<shared_ptr<Vertical>> & nodeList);

public:
    explicit DFD(std::filesystem::path const& path, char separator = ',', bool hasHeader = true);

    unsigned long long execute() override;
};
