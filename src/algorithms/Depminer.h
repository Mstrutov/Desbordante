#pragma once

#include "CSVParser.h"
#include "FDAlgorithm.h"
#include "depminer/util/CMAXGen.h"

class Depminer : public FDAlgorithm {
private:
    void lhsForColumn(const std::unique_ptr<Column>& column, std::vector<CMAXSet> const & cmaxSets);
    CMAXSet genFirstLevel(std::vector<CMAXSet> const & cmaxSets, Column attribute, std::unordered_set<Vertical>& li);
    std::unordered_set<Vertical> genNextLevel(std::unordered_set<Vertical> const& li);
    bool checkJoin(Vertical const& _p, Vertical const& _q);
public:
    
    explicit Depminer(std::filesystem::path const& path, char separator = ',', bool hasHeader = true) 
        : FDAlgorithm(path, separator, hasHeader) {}

    unsigned long long execute() override;
    std::unique_ptr<ColumnLayoutRelationData> relation;
    const RelationalSchema* schema;
};
