//
// Created by alexandrsmirn
//

#pragma once

#include "FDAlgorithm.h"

class DFD : public FDAlgorithm {
private:
    //TODO
public:
    explicit DFD(std::filesystem::path const& path, char separator = ',', bool hasHeader = true)
        : FDAlgorithm(path, separator, hasHeader) {}

    unsigned long long execute() override;
};
