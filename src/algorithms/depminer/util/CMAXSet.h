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

private:
    std::shared_ptr<Column> column;
    std::set<Vertical> columnCombinations;
public:
    CMAXSet(std::shared_ptr<Column> column) : column(column){};
    CMAXSet() = default;
    void addCombination(Vertical combination){columnCombinations.insert(combination);}
    std::set<Vertical> getCombinations(){return columnCombinations;}
    std::shared_ptr<Column> getColumn() const{ return column; }
    bool operator<(CMAXSet const& rhs) const{
        return *(this->getColumn()) < *(rhs.getColumn());
    }
};
