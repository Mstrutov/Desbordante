#pragma once

#include <unordered_set>

// #include "ColumnCombination.h"
// #include "ColumnData.h"
#include "ColumnLayoutRelationData.h"
// #include "RelationalSchema.h"
#include "custom/CustomHashes.h"

class CMAXSet{
private:
    Column column;
    std::unordered_set<Vertical> columnCombinations;
public:
    CMAXSet(Column column) : column(column) {};
    void makeNewCombinations(std::unordered_set<Vertical> comb){
        this->columnCombinations = comb;
    }
    void addCombination(Vertical const & combination){ columnCombinations.insert(combination); }
    const std::unordered_set<Vertical> & getCombinations() const { return columnCombinations; }
    const Column & getColumn() const { return column; }
};
