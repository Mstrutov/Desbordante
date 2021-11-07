#pragma once

#include <chrono>
#include <iostream>
#include <iomanip>
#include <list>
#include <memory>

#include "ColumnCombination.h"
#include "ColumnData.h"
#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"

class MAXSet{
private:
    Column column;
    std::set<Vertical> columnCombinations;
public:
    MAXSet(Column column) : column(column) {};
    MAXSet() = default;
    bool operator<(MAXSet const& rhs) const{
        return this->getColumn() < rhs.getColumn();
    }
    void makeNewCombinations(std::set<Vertical> comb){
        this->columnCombinations = comb;
    }
    void addCombination(Vertical combination){ columnCombinations.insert(combination); }
    std::set<Vertical> getCombinations(){ return columnCombinations; }
    Column getColumn() const{ return column; }
};
