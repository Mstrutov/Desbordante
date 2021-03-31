#pragma once

#include <filesystem>
#include <list>

#include "CSVParser.h"
#include "FD.h"

class FDAlgorithm {
protected:
    CSVParser inputGenerator_;
    std::list<FD> fdCollection_;
public:
    explicit FDAlgorithm (std::filesystem::path const& path, char separator = ',', bool hasHeader = true)
            : inputGenerator_(path, separator, hasHeader) {}

    virtual void registerFD(FD fdToRegister) { fdCollection_.push_back(std::move(fdToRegister)); }
    virtual void registerFD(Vertical lhs, Column rhs) { fdCollection_.emplace_back(std::move(lhs), std::move(rhs)); }
    virtual std::string getJsonFDs();
    virtual unsigned int fletcher16();
    virtual unsigned long long execute() = 0;
    virtual ~FDAlgorithm() = default;
};
