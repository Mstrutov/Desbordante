#pragma once

#include <random>
#include <stack>

#include "FDAlgorithm.h"
#include "Vertical.h"
#include "DFD/PartitionStorage/PartitionStorage.h"

class DFD : public FDAlgorithm {
private:
    std::unique_ptr<PartitionStorage> partitionStorage;
    std::unique_ptr<ColumnLayoutRelationData> relation;
    std::vector<Vertical> uniqueColumns;

    std::mutex mutable registerFdMutex;
    unsigned int numberOfThreads;

    void registerFD(Vertical vertical, Column rhs) override;

public:
    explicit DFD(std::filesystem::path const& path, char separator = ',', bool hasHeader = true, unsigned int parallelism = 0);

    unsigned long long execute() override;
};
