#pragma once

#include "CSVParser.h"
#include "FDAlgorithm.h"
#include "PliBasedFDAlgorithm.h"
#include "depminer/CMAXSet.h"

class Depminer : public PliBasedFDAlgorithm {
private:
    void lhsForColumn(std::unique_ptr<Column> const& column, std::vector<CMAXSet> const& cmaxSets);
    static CMAXSet genFirstLevel(std::vector<CMAXSet> const& cmaxSets, Column const& attribute, std::unordered_set<Vertical>& li);
    static std::unordered_set<Vertical> genNextLevel(std::unordered_set<Vertical> const& li);
    static bool checkJoin(Vertical const& _p, Vertical const& _q);
    std::vector<CMAXSet> generateCMAXSets(std::unordered_set<Vertical> const& agreeSets);
    double progressStep = 0;
public:
    
    explicit Depminer(std::filesystem::path const& path, char separator = ',', bool hasHeader = true) 
        : PliBasedFDAlgorithm(path, separator, hasHeader, true, {"AgreeSets generation", "Finding CMAXSets", "Finding LHS"}) {}

    unsigned long long executeInternal() override;
    const RelationalSchema* schema = nullptr;
};
