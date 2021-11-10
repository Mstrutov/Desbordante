#pragma once

#include "CSVParser.h"
#include "FDAlgorithm.h"
#include "depminer/util/CMAXGen.h"

class Depminer : public FDAlgorithm {
private:
    CMAXSet genFirstLevel(std::set<CMAXSet> cmaxSets, Column attribute, std::unordered_set<Vertical>& li);
    std::unordered_set<Vertical> genNextLevel(std::unordered_set<Vertical> const& li);
public:
    
    explicit Depminer(std::filesystem::path const& path, char separator = ',', bool hasHeader = true) 
        : FDAlgorithm(path, separator, hasHeader) {}

    unsigned long long execute() override;

};
