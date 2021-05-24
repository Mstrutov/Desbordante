#pragma once

#include "CSVParser.h"
#include "FDAlgorithm.h"
#include "depminer/util/CMAXSet.h"

class Depminer : public FDAlgorithm {
private:
    CMAXSet genFirstLi(std::set<CMAXSet> cmaxSets, Column attribute, std::set<Vertical>& li);
    std::set<Vertical> genNextLi(std::set<Vertical> const& li);
public:
    
    explicit Depminer(std::filesystem::path const& path, char separator = ',', bool hasHeader = true) 
        : FDAlgorithm(path, separator, hasHeader) {}

    unsigned long long execute() override;

};
