#pragma once

#include <unordered_set>

#include "ColumnLayoutRelationData.h"
#include "custom/CustomHashes.h"

class CMAXSet{
private:
    Column column;
    std::unordered_set<Vertical> columnCombinations;
public:
    CMAXSet(Column const& column) : column(column) {};
    void makeNewCombinations(std::unordered_set<Vertical> comb) {
        this->columnCombinations = comb;
    }
    void addCombination(Vertical const& combination) { columnCombinations.insert(combination); }
    std::unordered_set<Vertical> const& getCombinations() const { return columnCombinations; }
    Column const& getColumn() const { return column; }
};
