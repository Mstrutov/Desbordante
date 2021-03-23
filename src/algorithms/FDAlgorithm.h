#pragma once

#include "CSVParser.h"
#include <filesystem>

class FDAlgorithm {
protected:
    CSVParser inputGenerator_;
public:
    explicit FDAlgorithm (std::filesystem::path const& path, char separator = ',', bool hasHeader = true)
            : inputGenerator_(path, separator, hasHeader) {}
    virtual unsigned long long execute() = 0;
};