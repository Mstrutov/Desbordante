//
// Created by alexandrsmirn
//

#pragma once

#include <random>

#include "FDAlgorithm.h"

class DFD : public FDAlgorithm {
private:
    std::random_device rd; // для генерации случайных чисел
    std::mt19937 gen;

    void findLHSs(shared_ptr<Column> rhs, shared_ptr<RelationalSchema> schema); //TODO: нужен ли второй параметр?; мб переименовать типа findDeps
    shared_ptr<Vertical> randomWalk(std::list<shared_ptr<Vertical>> const& verticalList);
public:
    explicit DFD(std::filesystem::path const& path, char separator = ',', bool hasHeader = true)
        : FDAlgorithm(path, separator, hasHeader), gen(rd()) {}

    unsigned long long execute() override;
};
