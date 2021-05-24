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

class CMAXSet{

protected:
    Column column;
    std::set<Vertical> columnCombinations;
public:
    CMAXSet(Column column) : column(column){};
    void addCombination(Vertical combination){columnCombinations.insert(combination);}
    std::set<Vertical> getCombinations(){return columnCombinations;}
    Column getColumn() const{ return column; }
    bool operator<(CMAXSet const& rhs) const{
        return this->getColumn() < rhs.getColumn();
    }
};
