#pragma once

#include <chrono>
#include <iostream>
#include <iomanip>
#include <list>
#include <memory>
#include <unordered_set>

#include <boost/functional/hash.hpp>

#include "ColumnCombination.h"
#include "ColumnData.h"
#include "ColumnLayoutRelationData.h"
#include "RelationalSchema.h"
#include "custom/CustomHashes.h"

class MAXSet{
private:
    Column column;
    std::unordered_set<Vertical> columnCombinations;
public:
    MAXSet(Column column) : column(column) {};
    MAXSet() = default;
    bool operator<(MAXSet const& rhs) const{
        return this->getColumn() < rhs.getColumn();
    }
    void makeNewCombinations(std::unordered_set<Vertical> comb){
        this->columnCombinations = comb;
    }
    void addCombination(Vertical combination){ columnCombinations.insert(combination); }
    std::unordered_set<Vertical> getCombinations(){ return columnCombinations; }
    Column getColumn() const{ return column; }
};
