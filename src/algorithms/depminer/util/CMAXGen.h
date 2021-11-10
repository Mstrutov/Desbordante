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
#include "MAXSet.h"

using CMAXSet = MAXSet;

class CMAXGen{
private:
    std::vector<CMAXSet> cmaxSets;
    const RelationalSchema* schema;
public:
    CMAXGen(const RelationalSchema* schema) : schema(schema){};
    ~CMAXGen() = default;
    void execute(std::unordered_set<Vertical>& agreeSets);
    std::vector<CMAXSet> getCmaxSets(){
        return this->cmaxSets;
    }
};