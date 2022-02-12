#pragma once

#include <random>
#include <stack>

#include "PliBasedFDAlgorithm.h"
#include "Vertical.h"
#include "DFD/PartitionStorage/PartitionStorage.h"

class DFD : public PliBasedFDAlgorithm {
private:
    std::unique_ptr<PartitionStorage> partition_storage_;
    std::vector<Vertical> unique_columns_;

    unsigned int number_of_threads_;

    unsigned long long ExecuteInternal() override;

public:
    explicit DFD(std::filesystem::path const& path,
                 char separator = ',', bool has_header = true,
                 unsigned int parallelism = 0);
    explicit DFD(std::shared_ptr<ColumnLayoutRelationData> relation, unsigned int parallelism = 0);
};
