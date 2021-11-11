#pragma once

#include "CSVParser.h"
#include "FDAlgorithm.h"
#include "depminer/MAXSet.h"

class Depminer : public FDAlgorithm {
private:
    void lhsForColumn(const std::unique_ptr<Column>& column, std::vector<CMAXSet> const & cmaxSets);
    CMAXSet genFirstLevel(std::vector<CMAXSet> const & cmaxSets, Column attribute, std::unordered_set<Vertical>& li);
    std::unordered_set<Vertical> genNextLevel(std::unordered_set<Vertical> const& li);
    bool checkJoin(Vertical const& _p, Vertical const& _q);
    std::vector<CMAXSet> generateCMAXSets(std::unordered_set<Vertical>& agreeSets);
    double progressStep;
public:
    
    explicit Depminer(std::filesystem::path const& path, char separator = ',', bool hasHeader = true) 
        : FDAlgorithm(path, separator, hasHeader, {"AgreeSets generation", "Finding CMAXSets", "Finding LHS"}) {}

    unsigned long long execute() override;
    std::unique_ptr<ColumnLayoutRelationData> relation;
    const RelationalSchema* schema;
};
