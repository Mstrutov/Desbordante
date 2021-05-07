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
#include "CMAXSet.h"
#include "MAXSet.h"

class CMAXGen{
private:
    std::set<CMAXSet> cmaxSets;
    std::set<MAXSet> maxSets;
    std::shared_ptr<RelationalSchema> schema;
    void MaxSetsGenerate(std::set<Vertical> agreeSets);
    void CMaxSetsGenerate();
public:
    CMAXGen(std::shared_ptr<RelationalSchema> schema) : schema(schema){};
    CMAXGen() = default;
    ~CMAXGen() = default;
    void execute(std::set<Vertical> agreeSets);
};